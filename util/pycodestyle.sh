#!/usr/bin/env bash

pycodestyle --max-line-length=120 --verbose "${@}"
