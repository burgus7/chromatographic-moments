import pandas as pd
import math
import matplotlib.pyplot as plt

# add data to a dataframe
data = pd.read_table("data.csv", sep='\t')

# isolate concentration of B, conductivity, and time
concentration_B = data["percent_B"]
conductivity = data["mS/cm"]
time = data["min"]

# verify that the time intervals are approximately similar and store concentrations of B
# intervals stored as delta_t
all_conc_B = []
delta_t = []
for i in range(1, time.size):
    cur_time = time[i]
    prev_time = time[i-1]
    delta_t.append(cur_time - prev_time)
    if concentration_B[i] != concentration_B[i - 1]:
        all_conc_B.append(concentration_B[i])


# convert from list to dataframe, mainly because I'm learning Pandas right now and want more practice with it.
# thanks Han for suggesting that I learn Python when we met last semester :)
del_t = pd.Series(delta_t)
del_t_standard_error = del_t.std() / math.sqrt(del_t.size)
print(f"The 95% confidence interval for dt is: {round(del_t.mean(), 4)} +/- {round(del_t_standard_error * 1.96, 4)}.")
print(f"Since the standard error ({round(del_t_standard_error, 4)} is statistically insignificant, "
      f"the average delta t is used instead of finding the actual difference for each time step.")

dt = del_t.mean()

''' Moment Analysis Using Equations - Conductivity '''

salt_moments = []
start_index = 0
for salt_num in range(1, 5):
    # zeroth moment is integral(conductivity * dt), stored as M0
    M0 = 0
    for i in range(start_index, conductivity.size-1):
        cond = conductivity[i]
        rect_area = cond * dt
        M0 += rect_area
        if time[i + 1] > salt_num * 10:
            break
    M0 = round(M0, 4)

    # first moment is integral(conductivity * t * dt)/M0, stored as M1
    M1 = 0
    for i in range(start_index, conductivity.size-1):
        cond = conductivity[i]
        t = time[i]
        rect_area = cond * t * dt
        M1 += rect_area
        if time[i + 1] > salt_num * 10:
            break
    M1 = round(M1 / M0, 4)

    # second moment is integral(conductivity * (t - M1) ^ 2 * dt)/M0, stored as M2
    M2 = 0
    for i in range(start_index, conductivity.size-1):
        cond = conductivity[i]
        t = time[i]
        rect_area = cond * (t - M1) * (t - M1) * dt
        M2 += rect_area
        if time[i + 1] > salt_num * 10:
            start_index = i + 1
            break
    M2 = round(M2 / M0, 4)

    # store salt moments in list temporarily
    salt_moments.append([salt_num, f"{all_conc_B[salt_num-1]}%", M0, M1, M2])

moments = pd.DataFrame(salt_moments, columns=['Salt No.', 'Percent B', 'M0', 'M1', 'M2'])
print(moments)

''' Moment Analysis Using Transitions - Derivative of Conductivity '''

derivative_conductivity = []
for i in range(conductivity.size-1):
    delta_conductivity = conductivity[i+1] - conductivity[i]
    delta_time = del_t[i]
    slope = delta_conductivity/delta_time
    derivative_conductivity.append(slope)

# account for the off-by-one error caused by taking the difference
derivative_conductivity.append(slope)

# plot derivative
plt.plot(time, concentration_B, 'b--', label="% B")
plt.plot(time, conductivity, 'k-', label="Conductivity")
plt.plot(time, derivative_conductivity, 'go-', label="Derivative")
plt.vlines(0, 0, 0, linestyles ="dotted", colors ="red", label="First Moments")
# plot horizontal lines for the moments
for i in range(0, 4):
    x = moments.iloc[i]["M1"]
    plt.vlines(x, 0, 30, linestyles ="dotted", colors ="red")
plt.legend()
plt.title("Conductivity Profile")
plt.xlabel("Time (min)")
plt.ylabel("% Buffer B\n"
           "Conductivity (mS/cm)\n"
           "Derivative of Conductivity (mS/(cm * min))")
plt.grid()
plt.show()

