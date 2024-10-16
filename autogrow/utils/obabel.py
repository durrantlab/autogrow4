import os

from autogrow.utils.logging import log_warning


def obabel_convert_cmd(
    in_file: str, out_file: str, obabel_path: str, extra_params: str = ""
) -> str:
    # prot_and_3d: bool = False
    in_ext = in_file.split(".")[-1]
    out_ext = out_file.split(".")[-1]

    if in_ext == "vina":
        in_ext = "pdbqt"

    cmd = f'{obabel_path} -i{in_ext} "{in_file}" -o{out_ext}'
    if extra_params != "":
        # cmd += " --gen3d --p 7.4"
        cmd += f" {extra_params}"
    cmd += f' -e -O "{out_file}"'
    cmd += "  > /dev/null 2>&1"

    return cmd


def obabel_convert(
    in_file: str, out_file: str, obabel_path: str, extra_params: str = ""
) -> bool:
    cmd = obabel_convert_cmd(in_file, out_file, obabel_path, extra_params)
    try:
        os.system(cmd)
    except Exception as e:
        log_warning(f"Could not convert with obabel: {in_file}")
        return False

    if not os.path.exists(out_file):
        log_warning(f"Could not convert with obabel: {in_file}")
        return False

    with open(out_file, "r") as f:
        content = f.read().strip()
        if content == "":
            log_warning(f"Could not convert with obabel: {in_file}")
            return False

    return True
