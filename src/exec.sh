#!/bin/bash
echo "Start: $1"
nohup $1 > /dev/null 2>&1
echo "Finish: $1"