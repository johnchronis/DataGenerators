add_library(zipf zipf.h zipf.cpp)

add_executable(datagen datagen.cpp)

target_link_libraries(datagen zipf absl::flags absl::flags_parse)
