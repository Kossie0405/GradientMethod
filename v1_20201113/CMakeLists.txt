cmake_minimum_required(VERSION 2.8)

project(Main) # プロジェクト名を決める 
set(CMAKE_CXX_STANDARD 14) # C++14で書くよCMakeに教えてあげる 
include_directories(/usr/include/octave-5.2.0) #インクルードディレクトリ指定 
link_directories(/usr/lib/x86_64-linux-gnu) #ライブラリディレクトリ指定 
add_executable(optv2 optv2.cpp) #吐き出す実行ファイルと使うソースコードの指定 
target_link_libraries(optv2 octave) # 使うライブラリの指定。add_executableの後に書かないとダメ。
