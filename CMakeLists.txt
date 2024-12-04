cmake_minimum_required(VERSION 3.5)

project(MyQtApp VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Update this path based on where you found Qt5Config.cmake
set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/qt/6.7.2_1/lib/cmake")

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_executable(MyQtApp main.cpp)
target_link_libraries(MyQtApp Qt6::Widgets)

