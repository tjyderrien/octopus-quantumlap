#!/usr/bin/env python3
## Change character in .git/config (equivalent of sed)
# This comes in frame of change of repository names introduced on March 2nd 2021
# octopus-slabs --> octopus_slabs
# spp-extended-theory --> spp_extended_theory
# Contact: <derrien@fzu.cz>

Header = "RenameRepository.py: "

import shutil as su

def ReplaceInConfigFile(FormerName, NewName):  # {{{

    # Read the content of file
    try:
        f1 = open('.git/config', 'r')
        input_data = f1.read()
        f1.close()
    except:
        print(Header + "** Please execute this script at the root of the repository.)
        exit()

        # backup the config file in case
    try:
        su.copyfile('.git/config', '.git/config.20210302.bak')
    except:
        print(Header + "** Error in backup of the CONFIG file. STOP to prevent errors. Check the folder where you "
                       "execute this. ")
        exit()

    # Replace the target string
    input_data = input_data.replace( FormerName, NewName)

    # Write the output to same file
    f2 = open('.git/config', 'w')
    f2.write(input_data)
    f2.close()


# }}}

ReplaceInConfigFile('octopus-slabs', 'octopus_slabs')
ReplaceInConfigFile('spp-extended-theory', 'spp_extended_theory')