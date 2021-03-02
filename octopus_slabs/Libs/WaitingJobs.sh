#!/bin/bash
source machines.env
CheckComputerName
${QueueCommand} | grep 'tjyderri'
