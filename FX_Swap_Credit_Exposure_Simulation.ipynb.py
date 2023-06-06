#!/usr/bin/env python
# coding: utf-8

# # FX Swap Credit Exposure Simulation
# 
# In this notebook, we will simulate the credit exposures for three FX swap contracts. The counterparties are A, B, and C, and the notional amounts are 3M, 6M, and 8M US Dollars with JPY, INR and EUR for the tenor of 5,6,7 years respectively.
# 
# We will make the following assumptions:
# 
# 1. The FX rates and interest rates follow a geometric Brownian motion. This means that the rates will have a certain drift (average rate of change) and volatility (standard deviation of the rate of change).
# 2. We will assume some values for the drift and volatility of the FX rates and interest rates.
# 3. For simplicity, we will assume that the interest rates in all the currencies are the same.
# 4. We will also need to assume the recovery rate for the calculation of the Credit Value Adjustment (CVA).
# 5. For the calculation of the Potential Future Exposure (PFE), we will need to assume a confidence level.
# 6. We will assume that the swap contracts are plain vanilla swaps with semi-annual payments.
# 
# The steps are as follows:
# 
# 1. Simulate the FX rates and interest rates using the geometric Brownian motion.
# 2. Calculate the leg payments for each contract for each year.
# 3. Calculate the Expected Exposure (EE), Expected Positive Exposure (EPE), and Maximum Potential Future Exposure (MPFE) for each year.
# 4. Calculate the Credit Value Adjustment (CVA) for each contract.
# 5. Plot the exposures and leg payments for each contract.
# 6. Prepare a final report with all the results.

# In[ ]:


# Import necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

# Set the seed for reproducibility
np.random.seed(0)

# Define the parameters for the simulation
notional_amounts = [3e6, 6e6, 8e6]  # Notional amounts in USD
tenors = [5, 6, 7]  # Tenors in years
currencies = ['JPY', 'INR', 'EUR']  # Counterparty currencies
fx_rates = [110, 75, 0.85]  # Initial FX rates (USD per foreign currency)
interest_rate = 0.02  # Interest rate (assumed to be the same for all currencies)
volatility = 0.15  # Volatility of the FX rates
drift = 0.01  # Drift of the FX rates
recovery_rate = 0.4  # Recovery rate for the CVA calculation
confidence_level = 0.95  # Confidence level for the PFE calculation
num_simulations = 1000  # Number of simulations
num_intervals = 100  # Number of intervals per year for the simulation

# In[ ]:


# Define a function to simulate the FX rates
def simulate_fx_rates(fx_rate, drift, volatility, tenor, num_simulations, num_intervals):
    dt = tenor / num_intervals
    fx_rates = np.zeros((num_simulations, num_intervals + 1))
    fx_rates[:, 0] = fx_rate
    for t in range(1, num_intervals + 1):
        brownian = np.random.standard_normal(num_simulations)
        fx_rates[:, t] = fx_rates[:, t - 1] * np.exp((drift - 0.5 * volatility**2) * dt + volatility * np.sqrt(dt) * brownian)
    return fx_rates

# Define a function to calculate the leg payments
def calculate_leg_payments(notional_amount, fx_rates, interest_rate, tenor):
    num_payments = tenor * 2  # Semi-annual payments
    payments = np.zeros((fx_rates.shape[0], num_payments))
    for t in range(num_payments):
        payments[:, t] = notional_amount * fx_rates[:, int(t / 2 * fx_rates.shape[1] / num_payments)] * interest_rate / 2
    return payments

# Define a function to calculate the exposures
def calculate_exposures(leg_payments):
    exposures = np.cumsum(leg_payments, axis=1)
    return exposures

# Define a function to calculate the CVA
def calculate_cva(exposures, recovery_rate):
    cva = (1 - recovery_rate) * np.mean(np.max(exposures, axis=1))
    return cva

# Define a function to calculate the PFE
def calculate_pfe(exposures, confidence_level):
    pfe = np.percentile(exposures, confidence_level * 100, axis=0)
    return pfe

# In[ ]:


# Run the simulation for each contract
for i in range(3):
    print(f'Running simulation for contract with counterparty {chr(65 + i)}...')

    # Simulate the FX rates
    fx_rates_simulated = simulate_fx_rates(fx_rates[i], drift, volatility, tenors[i], num_simulations, num_intervals * tenors[i])

    # Calculate the leg payments
    leg_payments = calculate_leg_payments(notional_amounts[i], fx_rates_simulated, interest_rate, tenors[i])

    # Calculate the exposures
    exposures = calculate_exposures(leg_payments)

    # Calculate the CVA
    cva = calculate_cva(exposures, recovery_rate)
    print(f'The CVA for contract with counterparty {chr(65 + i)} is {cva}')

    # Calculate the PFE
    pfe = calculate_pfe(exposures, confidence_level)
    print(f'The PFE for contract with counterparty {chr(65 + i)} is {pfe}')

    # Plot the exposures
    plt.figure(figsize=(10, 6))
    plt.plot(np.mean(exposures, axis=0), label='Expected Exposure')
    plt.plot(np.max(exposures, axis=0), label='Maximum Exposure')
    plt.plot(pfe, label='Potential Future Exposure')
    plt.title(f'Exposures for contract with counterparty {chr(65 + i)}')
    plt.xlabel('Time (semi-annual periods)')
    plt.ylabel('Exposure (USD)')
    plt.legend()
    plt.show()

    # Plot the leg payments
    plt.figure(figsize=(10, 6))
    plt.plot(np.mean(leg_payments, axis=0), label='Expected Leg Payments')
    plt.plot(np.max(leg_payments, axis=0), label='Maximum Leg Payments')
    plt.title(f'Leg payments for contract with counterparty {chr(65 + i)}')
    plt.xlabel('Time (semi-annual periods)')
    plt.ylabel('Leg Payments (USD)')
    plt.legend()
    plt.show()

# In[ ]:


# Additional visualizations

# Plot the simulated FX rates for each contract
for i in range(3):
    fx_rates_simulated = simulate_fx_rates(fx_rates[i], drift, volatility, tenors[i], num_simulations, num_intervals * tenors[i])
    plt.figure(figsize=(10, 6))
    plt.plot(fx_rates_simulated.T)
    plt.title(f'Simulated FX rates for contract with counterparty {chr(65 + i)}')
    plt.xlabel('Time (semi-annual periods)')
    plt.ylabel('FX Rate (USD per {currencies[i]})')
    plt.show()

# Plot the distribution of the maximum exposure for each contract
for i in range(3):
    fx_rates_simulated = simulate_fx_rates(fx_rates[i], drift, volatility, tenors[i], num_simulations, num_intervals * tenors[i])
    leg_payments = calculate_leg_payments(notional_amounts[i], fx_rates_simulated, interest_rate, tenors[i])
    exposures = calculate_exposures(leg_payments)
    plt.figure(figsize=(10, 6))
    plt.hist(np.max(exposures, axis=1), bins=50, alpha=0.75)
    plt.title(f'Distribution of Maximum Exposure for contract with counterparty {chr(65 + i)}')
    plt.xlabel('Maximum Exposure (USD)')
    plt.ylabel('Frequency')
    plt.show()
