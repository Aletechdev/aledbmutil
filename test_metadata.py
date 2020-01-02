from metadata import get_condition_val_str


input_list = ['a', 'b', 'c']
assert(get_condition_val_str(input_list)=="b c")
assert(input_list==['a', 'b', 'c'])