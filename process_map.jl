using CSV

name = "./map/cornea_map_20by20_001.dat"
df = CSV.read(name, dataFrame; header=false)