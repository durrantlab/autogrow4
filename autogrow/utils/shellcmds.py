import os
import sys


def determine_bash_timeout_vs_gtimeout() -> str:
    """
    This function tests whether we should use the BASH command "timeout" (for linux)
     or the coreutils function "gtimeout" for MacOS which can be obtained
     through homebrew

    Returns:
    :returns: str timeout_option: A string either "timeout" or "gtimeout" describing
     whether the bash terminal is able to use the bash function timeout or gtimeout
    """

    if sys.platform.lower() in ["linux", "linux2"]:
        # Should be true and default installed in all Linux machines
        return "timeout"

    command = 'timeout 1 echo " "'
    # Running the os.system command for command will return 0,1, or 32512
    # 0 means that the timeout function works (most likely this is a linux os)
    # 32512 means that the timeout function DOES NOT Work (most likely this is MacOS)

    try:  # timeout or gtimeout
        timeout_result = os.system(f"g{command}")
    except Exception as e:
        raise Exception(
            "Something is very wrong. This OS may not be supported \
            by Autogrow or you may need to execute through Bash."
        ) from e
    if timeout_result == 0:
        return "gtimeout"
    print("gtimeout failed to run, we will check timeout")

    try:  # timeout or gtimeout
        timeout_result = os.system(command)
    except Exception as exc:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
            Autogrow or you may need to execute through Bash."
        ) from exc

    if timeout_result == 0:
        return "timeout"
    printout = (
        "Need to install GNU tools for Bash to work. \n"
        + "This is essential to use Bash Timeout function in Autogrow. \n"
    )
    printout += "\t This will require 1st installing homebrew. \n"
    printout += "\t\t Instructions found at: https://brew.sh/ \n"
    printout += "\t Once brew is installed, please run:"
    printout += " sudo brew install coreutils \n\n"
    print(printout)
    raise Exception(printout)
