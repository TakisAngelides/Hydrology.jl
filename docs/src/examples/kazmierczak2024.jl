using Hydrology

float_type = Float64
path = "/Users/taange001/Documents/Coding/Hydrology.jl/src/input/Kazmierczak2024/THWAITES2km_m3_HAB_toto.mat"
data_loading_function = load_Kazmierczak2024
c = PhysicalConstants{float_type}()
HM = HydrologyModel(path, data_loading_function, c)

update_q!(HM)
update_N!(HM)

visualize_field(HM.fields, :q)
visualize_field(HM.fields, :N)
