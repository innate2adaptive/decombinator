#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Convenience wrapper for running decombinator directly from source tree."""


from src.decombinator.io import cli_args
from src.decombinator.pipeline import run


if __name__ == "__main__":

    input = cli_args()
    run(cli_args=input)
