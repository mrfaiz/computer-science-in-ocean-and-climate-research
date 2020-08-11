# With "#" we can make comments

# First we want to define a function

def id(x):
    return x

print('value',id(3))
print('another tester', id(5))


# We want to define the sum_{k=start}^{end}f(k)
def sum_recursive(start, end, f):
    result = 0.0
    if(start< end):
        result = f(start)+sum_recursive(start+1,end,f)
    if(start==end):
        result = result + f(start)
    return result



def sum_iterative(start,end,f):
    result = 0.0
    if(start <=end):
        for k in range(start,end+1):
            print('k',k)
            result = result + f(k)
    if(start > end):
        return 'no sum to be calculated'
    return result

# 1+2+3+...+10=55
# The recursive version could have the problem, that we get an error:
# "maximum recursion depth exceeded" and this version is also in general slower
# Hence we programmed an iterative variant.
print(sum_recursive(1,100,id))
print(sum_iterative(1,10,id))

# 1+...+n=n(n+1)/2
import matplotlib.pyplot as plt
plt.plot([1, 2, 3, 4])
plt.ylabel('some numbers')
plt.show()