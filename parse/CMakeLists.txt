cmake_minimum_required(VERSION 3.9)

project(main)

add_definitions(-O3 -Wall -std=c99)

add_executable(main parser_demo.c parser.c)

