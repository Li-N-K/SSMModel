import numpy as np
import pandas as pd
import argparse
import json
import warnings
warnings.filterwarnings('ignore')


def Bfactor(theta, theta_fc, theta_w):
    """Beta factor to scale potential evapotranspiration (PET) based on soil moisture."""
    if theta >= theta_fc:
        return 1
    elif theta <= theta_w:
        return 0
    else:
        return (theta - theta_w) / (theta_fc - theta_w)


def calculate_Z(thickness):
    """Calculate gravitational potential height (z) at each layer interface."""
    total_column_height = sum(thickness)
    return [total_column_height - sum(thickness[:i]) for i in range(len(thickness) + 1)]


def SSMM_leapfrog(precipitation, pot_evaporation, layer_thickness, soil_props, lower_boundary='gravitational'):
    """
    Simulates multi-layer soil moisture using leapfrog integration.

    Args:
        precipitation: Series of daily precipitation [mm/day].
        pot_evaporation: Series of daily PET [mm/day].
        layer_thickness: List of thickness for each soil layer [mm].
        soil_props: Dictionary with keys: 'Qs', 'Field Capacity', 'Wilting Point',
                    'Ks(mm/h)', 'Psi(mm)', 'b'.
        lower_boundary: One of 'gravitational', 'ground_water', 'no_flow'.

    Returns:
        DataFrame of daily volumetric soil moisture per layer [%].
    """

    # add dummy last day for simulation of last actual day
    last_index = precipitation.index[-1] + pd.Timedelta(days=1)
    precipitation.loc[last_index] = 0
    pot_evaporation.loc[last_index] = 0

    num_layers = len(layer_thickness)
    num_steps = len(precipitation)
    time_step = 24
    hps = 24. / time_step

    # Soil parameters
    theta_s = soil_props['Qs']
    theta_fc = soil_props['Field Capacity']
    theta_w = soil_props['Wilting Point']
    K_s = soil_props['Ks(mm/h)']
    psi_s = soil_props['Psi(mm)']
    b = soil_props['b']

    grav_z = calculate_Z(layer_thickness)
    theta_initial = np.full(num_layers, theta_fc)

    theta = np.zeros((num_steps * time_step, num_layers))
    flux_top = np.zeros((num_steps * time_step, num_layers))
    flux_bottom = np.zeros((num_steps * time_step, num_layers))
    theta[0] = theta_initial

    for j in range(1, num_steps):
        b_factor = Bfactor(theta[(j - 1) * time_step, 0], theta_fc, theta_w)
        flux_top[(j - 1) * time_step + 1: j * time_step + 1, 0] = -(
            (precipitation.iloc[j] - pot_evaporation.iloc[j] * b_factor) / layer_thickness[0]
        ) / time_step

        for i in range(time_step * (j - 1) + 1, time_step * j + 1):
            for layer in range(num_layers):
                psi = psi_s * (theta[i - 1, layer] / theta_s) ** -b
                K = K_s * hps * (theta[i - 1, layer] / theta_s) ** (2 * b + 3)

                if layer == num_layers - 1:
                    if lower_boundary == 'gravitational':
                        psi_next = psi
                        K_next = K_s * hps
                        dz = layer_thickness[layer]
                    elif lower_boundary == 'ground_water':
                        psi_next = -psi
                        K_next = K_s * hps
                        dz = layer_thickness[layer]
                    elif lower_boundary == 'no_flow':
                        K_next = -K
                    else:
                        raise ValueError(f"Invalid lower_boundary: {lower_boundary}")

                    flux_bottom[i, layer] = (-((K + K_next) / 2) * ((psi + grav_z[layer]) - (psi_next + grav_z[layer + 1]))) / dz ** 2

                else:
                    psi_next = psi_s * (theta[i - 1, layer + 1] / theta_s) ** -b
                    K_next = K_s * hps * (theta[i - 1, layer + 1] / theta_s) ** (2 * b + 3)
                    dz = layer_thickness[layer]
                    dz_next = layer_thickness[layer + 1]
                    flux_bottom[i, layer] = (-((K + K_next) / 2) * ((psi + grav_z[layer]) - (psi_next + grav_z[layer + 1]))) / dz ** 2
                    flux_bottom[i, layer] = max((theta[i - 1, layer + 1] - theta[i - 1, layer]) / 2 / time_step, flux_bottom[i, layer])
                    flux_top[i, layer + 1] = (flux_bottom[i, layer] * dz ** 2) / dz_next ** 2

            for layer in range(num_layers):
                delta_theta = -(flux_top[i, layer] - flux_bottom[i, layer])
                if i == 1:
                    # Bootstrap step: theta[i - 2] would wrap around to the
                    # (still unset) last row of the array on the very first
                    # substep, so fall back to a forward step here instead of
                    # the leapfrog central difference.
                    theta[i, layer] = theta[i - 1, layer] + delta_theta
                elif i % 2 == 1:
                    theta[i, layer] = theta[i - 2, layer] + 2 * delta_theta
                else:
                    theta[i, layer] = theta[i - 1, layer] + delta_theta
                theta[i, layer] = min(max(theta[i, layer], theta_w), theta_s)

        if j % 100 == 0:
            for i in range(time_step * (j - 1) + 1, time_step * j + 1):
                for layer in range(num_layers):
                    delta_theta = -(flux_top[i, layer] - flux_bottom[i, layer])
                    theta[i, layer] = theta[i - 1, layer] + delta_theta
                    theta[i, layer] = min(max(theta[i, layer], theta_w), theta_s)

    column_names = [f'layer{i + 1}' for i in range(num_layers)]
    df5 = pd.DataFrame(theta[time_step - 1::time_step, :], columns=column_names)
    df5 = df5.set_index(precipitation.index)

    return df5.iloc[:-1]  # Removing the last dummy day


# ANSI escape codes for colored output
RED = '\033[91m'
GREEN = '\033[92m'
RESET = '\033[0m'


def main():
    parser = argparse.ArgumentParser(description="Run the SSMM_leapfrog soil moisture model.")
    parser.add_argument('--precip', type=str, required=True, help='Path to CSV file containing daily precipitation.')
    parser.add_argument('--pet', type=str, required=True, help='Path to CSV file containing daily potential evapotranspiration.')
    parser.add_argument('--soil', type=str, required=True, help='Path to JSON file containing soil properties.')
    parser.add_argument('--lower_boundary', type=str, default='gravitational',
                         choices=['gravitational', 'ground_water', 'no_flow'],
                         help='Lower boundary condition to use.')
    parser.add_argument('--output', type=str, default='soil_moisture_output.csv', help='Path to save output CSV.')

    args = parser.parse_args()

    # Load inputs
    precip = pd.read_csv(args.precip, index_col=0, parse_dates=True).squeeze()
    pet = pd.read_csv(args.pet, index_col=0, parse_dates=True).squeeze()
    layer_thickness = [300, 350, 400, 450, 500]
    with open(args.soil, 'r') as f:
        soil_props = json.load(f)

    # Run model
    print(f'{RED}Running the simulation...{RESET}')
    output = SSMM_leapfrog(precip, pet, layer_thickness, soil_props, args.lower_boundary)
    output['Total_Soil_Moisture'] = (output * layer_thickness).sum(axis=1)

    # Save results
    output.to_csv(args.output)

    print(f'{GREEN}Simulation complete! Results saved to {args.output}{RESET}')


if __name__ == '__main__':
    main()
