add_test( [==[Low Level Block Test]==] /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test/test_planner [==[Low Level Block Test]==]  )
set_tests_properties( [==[Low Level Block Test]==] PROPERTIES WORKING_DIRECTORY /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test)
add_test( [==[Low Level Block Test JSON]==] /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test/test_planner [==[Low Level Block Test JSON]==]  )
set_tests_properties( [==[Low Level Block Test JSON]==] PROPERTIES WORKING_DIRECTORY /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test)
add_test( [==[Basic Segment Test]==] /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test/test_planner [==[Basic Segment Test]==]  )
set_tests_properties( [==[Basic Segment Test]==] PROPERTIES WORKING_DIRECTORY /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test)
add_test( [==[Basic Planner Test]==] /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test/test_planner [==[Basic Planner Test]==]  )
set_tests_properties( [==[Basic Planner Test]==] PROPERTIES WORKING_DIRECTORY /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test)
add_test( [==[Large Small Planner Test]==] /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test/test_planner [==[Large Small Planner Test]==]  )
set_tests_properties( [==[Large Small Planner Test]==] PROPERTIES WORKING_DIRECTORY /Users/eric/Documents/proj/trajectory/src/cmake-cli-test/test)
set( test_planner_TESTS [==[Low Level Block Test]==] [==[Low Level Block Test JSON]==] [==[Basic Segment Test]==] [==[Basic Planner Test]==] [==[Large Small Planner Test]==])
