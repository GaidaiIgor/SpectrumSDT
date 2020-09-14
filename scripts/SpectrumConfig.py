import os.path as path
from typing import List

class SpectrumConfig:
    """ Provides access to information in config file """
    def __init__(self, config_path):
        self.config_path = config_path
        with open(config_path, "r") as config_file:
            config_lines = config_file.readlines()
        self.params = {}
        for line in config_lines:
            line = line.split("!")[0].strip()  # remove comments
            if len(line) == 0:
                continue
            tokens = line.split("=")
            if len(tokens) != 2:
                raise Exception('Error in config format')
            key = tokens[0].strip()
            value = tokens[1].strip()
            self.params[key] = value

    def get_launch_mode(self) -> str:
        return self.params["mode"]

    def get_rovib_coupling(self) -> int:
        return int(self.params["rovib_coupling"])

    def get_fix_basis_jk(self) -> int:
        if "fix_basis_jk" in self.params:
            return int(self.params["fix_basis_jk"])
        else:
            return 0

    def get_j(self) -> int:
        return int(self.params["J"])

    def get_ks(self) -> List[int]:
        K_str = self.params["K"]
        if K_str == "all":
            J = self.get_j()
            k_start = (J + self.get_parity()) % 2
            return [k_start, J]

        K_str_tokens = K_str.split("..")
        if len(K_str_tokens) == 1:
            return [int(K_str_tokens[0]), int(K_str_tokens[0])]
        else:
            return [int(K_str_tokens[0]), int(K_str_tokens[1])]

    def get_parity(self) -> int:
        return int(self.params["parity"])

    def get_symmetry(self) -> int:
        return int(self.params["symmetry"])

    def get_basis_root_path(self) -> str:
        if "fix_basis_jk" in self.params:
            return self.params["basis_root_path"]
        else:
            return "-1"

    def get_basis_j(self) -> int:
        if self.get_fix_basis_jk() == 1:
            return int(self.params["basis_J"])
        else:
            return -1

    def get_basis_k(self) -> int:
        if self.get_fix_basis_jk() == 1:
            return int(self.params["basis_K"])
        else:
            return -1

    def get_solver(self) -> str:
        if "solver" in self.params:
            return self.params["solver"]
        else:
            return None

    def get_number_of_states(self) -> int:
        return int(self.params["num_states"])

    def get_ncv(self) -> int:
        if "ncv" in self.params:
            return int(self.params["ncv"])
        else:
            return int(self.get_number_of_states() * 2) # default value

    def get_grid_folder_path(self) -> str:
        grid_folder_path = self.params["grid_path"]
        return self.get_full_path(grid_folder_path)

    def get_root_path(self) -> str:
        return self.params["root_path"]

    def get_test_mode(self) -> str:
        if "test_mode" in self.params:
            return self.params["test_mode"]
        else:
            return ""

    def get_number_of_wfs_to_print(self) -> int:
        return 0

    def get_full_path(self, rel_path):
        """ Resolves given relative path with respect to path to the current config instance file"""
        return path.normpath(path.join(path.split(self.config_path)[0], rel_path))

