import domain as dom
import yaml


class Parameters:
    def __init__(self, parameter_file_path: str):
        with open(parameter_file_path, "r") as file:
            file_dict = yaml.safe_load(file)

        self.derivative_order = file_dict["derivatives"]["order"]
        if (self.derivative_order % 2 != 0):
            raise RuntimeError(f"Unsoported derivative order {
                               self.derivative_order}. Derivative orders must be divisible by 2")

        self.last_iter = file_dict["time"]["last-iter"]
        self.courant_factor = file_dict["time"]["courant-factor"]

        self.domain = dom.Domain(
            file_dict["domain"]["start"],
            file_dict["domain"]["end"],
            file_dict["domain"]["points"],
            int(self.derivative_order / 2),
            self.courant_factor
        )

        self.id_type = file_dict["initial-data"]["type"]

        if (self.id_type == "standing-wave"):
            self.standing_wave_A = file_dict["initial-data"]["standing-wave"]["A"]
            self.standing_wave_Kx = file_dict["initial-data"]["standing-wave"]["Kx"]
            self.standing_wave_Ky = file_dict["initial-data"]["standing-wave"]["Ky"]
            self.standing_wave_Kz = file_dict["initial-data"]["standing-wave"]["Kz"]
        else:
            self.standing_wave_A = None
            self.standing_wave_Kx = None
            self.standing_wave_Ky = None
            self.standing_wave_Kz = None

        self.c1 = file_dict["RKAB"]["c1"]
        self.c2 = file_dict["RKAB"]["c2"]
        self.c3 = file_dict["RKAB"]["c3"]
        self.c4 = file_dict["RKAB"]["c4"]
        self.c5 = file_dict["RKAB"]["c5"]
        self.c6 = file_dict["RKAB"]["c6"]
        self.c7 = file_dict["RKAB"]["c7"]
        self.c8 = file_dict["RKAB"]["c8"]
