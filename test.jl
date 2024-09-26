using Printf
using Infiltrator

data = [
    [2.333333e-04, 2.333333e-04, 3.000000e-06, 1.200000e+00],
    [2.713249e-04, 3.260851e-04, 3.000000e-06, 1.200000e+00],
    [2.439990e-04, 2.660970e-04, 3.000000e-06, 1.200000e+00],
    [2.153969e-04, 4.494829e-04, 3.000000e-06, 1.200000e+00]
]

# Writing to the .dat file
open("output.dat", "w") do file
    for row in data
        formatted_row = [@sprintf("%e", x) for x in row]
        println(file, join(formatted_row, ","))
        @infiltrate
    end
end