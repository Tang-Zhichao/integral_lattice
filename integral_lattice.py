import numpy as np
import collections
import itertools
from fractions import Fraction

class intergral_lattice:
    def __init__(self,init_form):
        self.intersection_form = init_form
        self.dim = self.__dim()
        self.disc = self.__disc()
        self.dual_group= self.__dual_group()
    
    def __dim(self):
        return len(self.intersection_form[0])

    def __disc(self):
        return int(round(np.linalg.det(self.intersection_form)))

    def __dual_group(self):
        potential_range = itertools.product(range(self.disc),repeat=self.dim)
        result = collections.deque()
        for potential_elment in potential_range:
            potential_elment_array = np.asarray(potential_elment)
            test = np.dot(potential_elment_array,self.intersection_form) % self.disc
            if np.all(test == 0):
                result.append(potential_elment_array)
        result_array = np.reshape(np.asarray(result),(-1,self.dim))
        return result_array

    def __factor(self,num):
        factor  = []
        for i in range(2 , num + 1):
            while (num % i == 0):
                factor.append(i)
                num = num // i
        return factor

    def __prime_list(self,num=-1):
        if (num == -1):
            num = self.disc
        return self.__factor(num)

    def __is_subarray(array,target):
        for sub_array in array:
            if (np.array_equal(sub_array,target) == True):
                return True
        return False

    def __spanned_by(self,arr,target):
        if (len(arr) != 0  and len(target) != 0):
            for i in itertools.product(range(self.disc),repeat=len(target)):
                arr_test = np.zeros(self.dim)
                for j in range(len(target)):
                    arr_test = arr_test + i[j] * target[j]
                arr_test = arr_test % self.disc
                if np.array_equiv(arr_test,arr):
                    return True
        return False

    def __array_replace_with(self,arr,target_0,target_1):
        result = collections.deque()
        for sub_arr in arr:
            if (np.array_equiv(sub_arr,target_0) == False):
                result.append(sub_arr)
            else:
                result.append(target_1)
        return result

    def __order_of(self,element):
        target = np.asarray(element)
        if (intergral_lattice.__is_subarray(self.dual_group,target) == False):
            return 0
        for i in range(1 , self.disc + 1):
            target_i = (i * target) % (self.disc)
            if np.all(target_i == 0):
               return i

    def __product_of_order(self,array):
        product = 1
        for sub_arr in array:
            product = (product) * (self.__order_of(element=sub_arr))
        return product
    
    def __frac_output_of_generators(self,array):
        result= collections.deque()
        for r in array:
            temp = collections.deque()
            for m in range(self.dim):
                temp.append(Fraction(r[m],self.disc))
            result.append(np.asarray(temp))
        result_array = np.reshape(np.asarray(result),(-1,self.dim))
        print(f"There are {len(result_array)} generators:")
        for l in result_array:
            print("[",end="")
            for j in range(len(l)):
                if (j != 0):
                    print(" ",end="")
                if (l[j].is_integer() == True):
                    print("{}".format(l[j].numerator),end="")
                else:
                    print("{}/{}".format(l[j].numerator,l[j].denominator),end="")
                if (j != len(l)-1):
                    print(" ",end="")
            print("]")
        print("with orders:")
        print(self.order_list(array))

    def __max_order_of(self,list):
        max_order = 0
        for index in range(len(list)):
            max_order = max(max_order,self.__order_of(element=list[index]))
        return max_order

    def __quot_with(self,target):
        expected_size = int(self.disc / self.__product_of_order(target))
        result = collections.deque()
        for arr in self.dual_group:
            arr_order = self.__order_of(arr)
            if (arr_order == 1):
                result.append(arr)
            elif (arr_order <= expected_size):
                if (self.__spanned_by(arr,target) == False):
                    if (all(self.__spanned_by(arr - term, target) == False) for term in result):
                        result.append(arr)
                    if (len(result) == expected_size):
                        break
        for i in range(len(result)):
            for arr in self.dual_group:
                order_of_i = self.__order_of(result[i])
                if (np.all(arr == 0)):
                    continue
                if (self.__order_of(arr) < order_of_i):
                    if (self.__spanned_by(arr - result[i],target) == True):
                        result = self.__array_replace_with(result,result[i],arr)
        return result
    
    def order_list(self,list):
        order_list = []
        for index in range(len(list)):
            order_list.append(self.__order_of(element=list[index]))
        return order_list
    
    def dual_group_generator(self):
        remain_list = collections.deque(self.dual_group)
        genrators = collections.deque()
        while (len(remain_list) != 1):
            max_order = self.__max_order_of(remain_list)
            for item in remain_list:
                if (self.__order_of(item) == max_order):
                    genrators.append(item)
                    remain_list = self.__quot_with(genrators)
                    break
        generators_array = np.reshape(np.asarray(genrators),(-1,self.dim))
        return generators_array
    
    def dual_group_generator_frac(self):
        self.__frac_output_of_generators(self.dual_group_generator())
