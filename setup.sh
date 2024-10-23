#!/bin/bash

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Function to display help
display_help() {
    echo -e "${YELLOW}Usage: bash setup.sh [Package_manager]${NC}"
    echo -e "   Package_manager    Using conda or mamba to install."
    echo -e "   -h                 Display this help message."
}

# Check for help option
if [[ "$1" == "-h" ]]; then
    display_help
    exit 0
fi

# Check for package manager argument
if [[ "$1" != "conda" && "$1" != "mamba" ]]; then
    echo -e "${RED}Error: Invalid input. Please specify 'conda' or 'mamba'.${NC}"
    display_help
    exit 1
fi

PACKAGE_MANAGER=$1

# 1. Check the operating system version
OS="$(uname)"
if [[ "$OS" != "Linux" ]]; then
    echo -e "${RED}Error: This script only supports Linux systems.${NC}"
    exit 1
fi

# 2. Output computer configuration information
echo -e "${YELLOW}Computer Configuration Information:${NC}"
echo "-------------------"
echo -e "User: ${GREEN}$(whoami)${NC}"
echo -e "Current Directory: ${GREEN}$(pwd)${NC}"

# Get CPU information
CPU_MODEL=$(lscpu | grep 'Model name' | awk -F ': ' '{print $2}' | xargs)
CPU_COUNT=$(lscpu | grep '^CPU(s):' | awk '{print $2}')
CORE_COUNT=$(lscpu | grep '^Core(s) per socket:' | awk '{print $4}')
THREAD_COUNT=$(lscpu | grep '^Thread(s) per core:' | awk '{print $4}')
CPU_FREQ=$(lscpu | grep 'CPU MHz' | awk -F ': ' '{print $2}' | xargs)

echo -e "CPU Model: ${GREEN}$CPU_MODEL${NC}"
echo -e "CPU(s): ${GREEN}$CPU_COUNT${NC}"
echo -e "Core(s) per socket: ${GREEN}$CORE_COUNT${NC}"
echo -e "^Thread(s) per core: ${GREEN}$THREAD_COUNT${NC}"
echo -e "CPU Frequency: ${GREEN}$CPU_FREQ MHz${NC}"
echo -e "Memory: ${GREEN}$(free -h | grep '^Mem:' | awk '{print $2}')${NC}"

echo -e "          "

# Check GPU information
if command -v nvidia-smi &> /dev/null; then
    echo -e "GPU Information:"
    GPU_INFO=$(nvidia-smi --query-gpu=name,memory.total,memory.free,memory.used,driver_version --format=csv,noheader)
    GPU_COUNT=$(echo "$GPU_INFO" | wc -l)

    if [[ "$GPU_COUNT" -eq 0 ]]; then
        echo -e "${RED}Error: No NVIDIA GPUs found.${NC}"
        exit 1
    fi

    echo -e "Number of GPUs: ${GREEN}$GPU_COUNT${NC}"
    echo "$GPU_INFO" | while IFS=',' read -r name mem_total mem_free mem_used driver_version; do
        echo -e "GPU Model: ${GREEN}$(echo "$name" | xargs)${NC}"
        echo -e "  Total Memory: ${GREEN}$(echo "$mem_total" | xargs)${NC}"
        echo -e "  Free Memory: ${GREEN}$(echo "$mem_free" | xargs)${NC}"
        echo -e "  Used Memory: ${GREEN}$(echo "$mem_used" | xargs)${NC}"
        echo -e "  Driver Version: ${GREEN}$(echo "$driver_version" | xargs)${NC}"
    done
    echo "-------------------"
else
    # Fallback to lspci if nvidia-smi is not available
    echo -e "${YELLOW}Checking for NVIDIA GPU using lspci...${NC}"
    if lspci | grep -i nvidia &> /dev/null; then
        echo -e "${RED}Error: NVIDIA driver not found. Please ensure you have an NVIDIA GPU and drivers installed.${NC}"
        exit 1
    else
        echo -e "${RED}Error: No NVIDIA GPUs detected.${NC}"
        exit 1
    fi
fi


sleep 3

# 3. Execute package manager command
if [[ ! -f environment.yml ]]; then
    echo -e "${RED}Error: environment.yml file not found.${NC}"
    exit 1
fi

echo -e "${YELLOW}Creating environment using $PACKAGE_MANAGER...${NC}"
$PACKAGE_MANAGER env create -f environment.yml

# 4. Execute pip install command
# activate space
$PACKAGE_MANAGER init
$PACKAGE_MANAGER activate space

# check environment name
CURRENT_ENV=$(basename "$CONDA_PREFIX")

if [ "$CURRENT_ENV" != "space" ]; then
    echo "Error: Current environment is not 'space'."
    exit 1
else
    echo "Successfully activated environment: $CURRENT_ENV"
fi

echo -e "${YELLOW}Installing Python packages...${NC}"
pip install --no-deps bokeh==3.4.2 stlearn==0.4.12

echo -e "${GREEN}Script execution completed.${NC}"