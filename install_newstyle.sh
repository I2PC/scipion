#!/bin/sh

# Newstyle Xmipp installation script
# Description: This script is intended to be the one which guides the installation process with the new buildout system.
# Author: Ignacio Foche Perez (ifoche@cnb.csic.es)

# Interface variables
TITLE="Xmipp binary generation script"
BLACK='\033[30m'
WHITE='\033[37m'
YELLOW='\033[33m'
RED='\033[31m'
BLUE='\033[36m'


#Some flags variables



welcome message(){
  echo -e "${BLACK}0000000000000000000000000000000000000000000000000001"              
  echo -e "${BLACK}0000000000000000P!!00000!!!!!00000!!0000000000000001"
  echo -e "${BLACK}000000000000P'  ${RED}.:==.           ,=;:.  ${BLACK}\"400000000001"
  echo -e "${BLACK}0000000000!  ${RED}.;=::.::,         .=:.-:=,.  ${BLACK}!000000001"
  echo -e "${BLACK}0000000P' ${RED}.=:-......::=      ::;.......-:=. ${BLACK}\"0000001"
  echo -e "${BLACK}0000000,${RED}.==-.........::;.   .;::.-.......-=:${BLACK}.a#00001"
  echo -e "${BLACK}0000000'${RED}.=;......--:.:.:=:.==::.:--:.:...,=- ${BLACK}!000001"
  echo -e "${BLACK}000000    ${RED}-=...:.:.:-:::::=;::.::.:.:..-:=-    ${BLACK}00001"
  echo -e "${BLACK}0000P       ${RED}==:.:-::. ${YELLOW}.aa.${RED} :::::::::.::=:      ${BLACK}\"4001"
  echo -e "${BLACK}00001        ${RED}:;:::::..${YELLOW}:#0:${RED} =::::::::::::        ${BLACK}j#01"
  echo -e "${BLACK}0001   ${YELLOW}aa _aas  _aa_  .aa, _a__aa_.  .a__aas,    ${BLACK}j#1"
  echo -e "${BLACK}"'0001   '"${YELLOW}"'4WU*!4#gW9!##i .##; 3##P!9#Ai .##U!!Q#_   '"${BLACK}"'j#1'
  echo -e "${BLACK}"'001    '"${YELLOW}"'3O.  :xO.  ]Oi .XX: ]O( .  X2 .XC;  :xX.   '"${BLACK}"'01'
  echo -e "${BLACK}"'001    '"${YELLOW}"'dU.  :jU.  %Ui .WW: ]UL,..aXf .Ums  jd*    '"${BLACK}"'01'
  echo -e "${BLACK}"'0001   '"${YELLOW}"'4W.  :dW. .%W1 :WW: %WVXNWO~  .#U*#WV!    _'"${BLACK}"'01'
  echo -e "${BLACK}"'0001        '"${RED}"'.............. '"${YELLOW}"'%#1  '"${RED}"'.... '"${YELLOW}"'.#A)        '"${BLACK}"'j01'
  echo -e "${BLACK}"'00001      '"${RED}"':=::-:::::;;;;: '"${YELLOW}"'301 '"${RED}"'::::: '"${YELLOW}"'.0A)       '"${BLACK}"'j#01'
  echo -e "${BLACK}0000L    ${RED}.;::.--::::::::;: ,,..:::::... .::.   ${BLACK}_d001"
  echo -e "${BLACK}000000  ${RED}:;:...-.-.:.:::=;   =;:::---.:....::,  ${BLACK}00001"
  echo -e "${BLACK}000000!${RED}:::.....-.:.::::-     :=:-.:.-......:= ${BLACK}!00001"
  echo -e "${BLACK}000000a  ${RED}=;.........:;        .:;........,;;  ${BLACK}a00001"
  echo -e "${BLACK}"'0000000La '"${RED}"'--:_....:=-           :=:...:_:--'"${BLACK}"'_a#000001'
  echo -e "${BLACK}"'0000000000a  '"${RED}"'-=;,:=              -;::=:-  '"${BLACK}"'a000000001'
  echo -e "${BLACK}"'00000000000Laaa '"${RED}"'-\aa              ar- '"${BLACK}"'aaa00000000001'
  echo -e "${BLACK}00000000000000#00000000aaaaaad000000#00#0#0000000001"
  echo ""
  echo -e "${WHITE}Welcome to $TITLE ${WHITE}"
  echo ""
  return 0;
}

#function that prints the script usage help
helpMessage(){
  echo -e ""
  welcomeMessage
  echo -e "###################################"
  echo -e '# XMIPP INSTALLATION SCRIPT USAGE #'
  echo -e "###################################"
  echo -e "${RED}NAME"
  echo -e "${WHITE}  install_newstyle.sh - XMIPP installation new generation script. "
  echo -e ""
  echo -e "${RED}SYNOPSIS"
  echo -e "${WHITE}  ./xmipp_binary_generator.sh ${BLUE}[OPTIONS]${WHITE}"
  echo -e ""
  echo -e "${RED}DESCRIPTION"
  echo -e "${WHITE}  Script that automates the XMIPP compilation process. When this script is executed, the compilation sequence starts, depending on the selected options. If no option is given, the script follows the sequence:"
  echo -e "1- untar the external libraries."
  echo -e "2- Compile one by one every external library"
  echo -e "3- Compile a little python"
  echo -e "4- Create virtual enviroment"
  echo -e "5- Install on this virtualenv the needed python modules"
  echo -e "6- Launch SConscript to compile Xmipp"
  echo -e "7- Run the unitary tests"
  echo -e ""
  echo -e "No option is mandatory. The default behaviour (without any given option) Following options are accepted:"
  echo -e ""



  echo -e "${BLUE}--platform=${YELLOW}<PLATFORM>${WHITE},${BLUE} -p ${YELLOW}<PLATFORM>${WHITE}"
  echo -e "    Selects the platform for which the binary is going to be built. Default behaviour is OpenSuSE"
  echo -e "${BLUE}--list-platforms${WHITE},${BLUE} -lp${WHITE}"
  echo -e "    Lists the possible platforms Xmipp3.0 can be built for"
  echo -e "${BLUE}--architecture=${YELLOW}<ARCH>${WHITE},${BLUE} -a ${YELLOW}<ARCH>${WHITE}"
  echo -e "    Selects the architecture for which the binary is going to be built. Default behaviour is i386"
  echo -e "${BLUE}--list-architectures${WHITE},${BLUE} -la${WHITE}"
  echo -e "    Lists the possible architectures Xmipp3.0 can be built for"
  echo -e "${BLUE}--branch=${YELLOW}<BRANCH>${WHITE},${BLUE} -b ${YELLOW}<BRANCH>${WHITE}"
  echo -e "    Decides which branch is going to be built from the Xmipp repository. Default behaviour is 3.0"
  echo -e "${BLUE}--remote-host=${YELLOW}<HOST_IP>${WHITE},${BLUE} -r ${YELLOW}<HOST_IP>${WHITE}"
  echo -e "    Sets the IP for the remote host that has the VM and will do the binary generation itself. Default behaviour is localhost"
  echo -e "${BLUE}--iface=${YELLOW}<IFACE>${WHITE},${BLUE} -i ${YELLOW}<IFACE>${WHITE}"
  echo -e "    Interface to select the options needed in the script to work. Default behaviour is a dialog-based interface"
  echo -e "${BLUE}--remote-user=${YELLOW}<USER>${WHITE},${BLUE} -u ${YELLOW}<USER>${WHITE}"
  echo -e "    Sets the user for the connection with the remote host. Default behaviour is the same user that runs the script" 
  echo -e "${BLUE}--list-ifaces${WHITE},${BLUE} -li${WHITE}"
  echo -e "    Lists the possible interfaces this script can be executed with"
  echo -e "${BLUE}--date-path${WHITE},${BLUE} -d${WHITE}"
  echo -e "    Decides path where a file containing the binaries generation date will be placed to be imported from Xmipp wiki"
  echo -e "${BLUE}--wakeup-time=${YELLOW}<WAKEUP_TIME>${WHITE},${BLUE} -w ${YELLOW}<WAKEUP_TIME>${WHITE}"
  echo -e "    Provides a number of CPUs for the compilation process. Default value is 2"
  echo -e "${BLUE}-j ${YELLOW}<NUM_CPUS>${WHITE}"
  echo -e "    Provides a number of CPUs for the compilation process. Default value is 2"
  echo -e "${BLUE}--help,${BLUE} -h${WHITE}"
  echo -e "    Shows this help message"
  echo -e ""
  echo -e "${YELLOW}<PLATFORM>"
  echo -e "${WHITE}  Execute this script with ${BLUE}-lp${WHITE} or ${BLUE}--list-platforms${WHITE} argument to know the different options you have for this field"
  echo -e ""
  echo -e "${YELLOW}<ARCH>"
  echo -e "${WHITE}  Execute this script with ${BLUE}-la${WHITE} or ${BLUE}--list-architectures${WHITE} argument to know the different options you have for this field"
  echo -e ""
  echo -e "${YELLOW}<HOST_IP>"
  echo -e "${WHITE}  The IP for connecting the machine where Virtual Machines are placed and binaries are going to be built (slave binary generation script must be in its PATH environment var). Default behaviour is localhost"
  echo -e ""
  echo -e "${YELLOW}<USER>"
  echo -e "${WHITE}  The user for connecting the machine where Virtual Machines are placed and binaries are going to be built (slave binary generation script must be in its PATH environment var). Default behaviour is the same user that runs the script"
  echo -e ""
  echo -e "${YELLOW}<BRANCH>"
  echo -e "${WHITE}  The branch, in the git repository where Xmipp is hosted, which is going to be built to produce the binaries"
  echo -e ""
  echo -e "${YELLOW}<WAKEUP_TIME>"
  echo -e "${WHITE}  Time in seconds the Virtual Machine will wake for the wakeup process"
  echo -e ""
  return 0;
}

