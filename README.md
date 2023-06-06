# FX-Swap-Credit-Exposure-Simulation

In this notebook, we will simulate the credit exposures for three FX swap contracts. The
counterparties are A, B, and C, and the notional amounts are 3M, 6M, and 8M US Dollars
with JPY, INR and EUR for the tenor of 5,6,7 years respectively.
# We will make the following assumptions:
1. The FX rates and interest rates follow a geometric Brownian motion. This means that
the rates will have a certain drift (average rate of change) and volatility (standard
deviation of the rate of change).
2. We will assume some values for the drift and volatility of the FX rates and interest rates.
3. For simplicity, we will assume that the interest rates in all the currencies are the same.
4. We will also need to assume the recovery rate for the calculation of the Credit Value
Adjustment (CVA).
5. For the calculation of the Potential Future Exposure (PFE), we will need to assume a
confidence level.
6. We will assume that the swap contracts are plain vanilla swaps with semi-annual
payments.

# The steps are as follows:
1. Simulate the FX rates and interest rates using the geometric Brownian motion.
2. Calculate the leg payments for each contract for each year.
3. Calculate the Expected Exposure (EE), Expected Positive Exposure (EPE), and
Maximum Potential Future Exposure (MPFE) for each year.
4. Calculate the Credit Value Adjustment (CVA) for each contract.
5. Plot the exposures and leg payments for each contract.
6. Prepare a final report with all the results.
