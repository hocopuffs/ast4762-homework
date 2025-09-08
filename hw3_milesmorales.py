#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Shiloh
#08 Sep 2025
#Homework 3
print("Problem 2")

import numpy as np
from hw3_support_functions_file import square

print("\nh")
test_square_1 = np.arange(10)
print(square(test_square_1))
#used np.arange to generate an array from 0 to 9

print("\ni")
test_square_2 = np.linspace(0, 24, 25).reshape(5,5)
print(square(test_square_2))
#wording of question was slightly ambiguous
#I assumed the result was supposed to be whole numbers instead of 25 evenly spaced decimals between 0 and 25

from hw3_support_functions_file import squareplot
#Had to refresh Kernel to get it to work
print("Problem 3")

squareplot(1,7,5, saveplot="plot_hw3_milesmorales.pdf")


# In[ ]:




