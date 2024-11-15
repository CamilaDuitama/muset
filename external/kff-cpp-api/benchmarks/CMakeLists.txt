
# add the executable
add_executable(reader_bench
    reader.cpp
)

target_link_libraries(reader_bench
    kff
)

