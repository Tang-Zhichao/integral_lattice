#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <Eigen>

template<int N>
class INTEGRAL_LATTICE{
    private:
    Eigen::Matrix<int,N,N> intersection_form;
    int dim;
    int disc;

    class GROUP_ELEMENT{
        public:
        Eigen::Matrix<int,1,N> coordinates;
        int order;
        GROUP_ELEMENT(Eigen::Matrix<int,1,N> input_vector, int input_order){
            coordinates = input_vector;
            order = input_order;
        }
        GROUP_ELEMENT(std::vector<int> input_vector, int input_order){
            int iuput_vector_in_array[N];
            std::copy(input_vector.begin(),input_vector.end(),iuput_vector_in_array);
            coordinates = Eigen::Map<Eigen::Matrix<int,1,N>>(iuput_vector_in_array);
            order = input_order;
        }
        bool operator<(const GROUP_ELEMENT& other){
            if (this->order < other.order) return true;
            else return false;
        }
    };

    std::deque<GROUP_ELEMENT> dual_group;
    std::deque<GROUP_ELEMENT> dual_group_generators;
    
    template<typename T>
    int the_order_of(const T& input_vector){
        for (int i = 1; i <= disc; i++){
            int t = 0;
            for (int k = 0; k < N; k++)
                if ((input_vector[k] * i) % disc == 0) t++;
            if (t == N) return i;
        }
        return 0;
    }
    
    template<typename T>
    bool integral_check(const T& input_vector){
        for(int i = 0; i < N; i++){
            int t = 0;
            for(int j = 0; j < N; j++) t += input_vector[j] * intersection_form(i,j);
            if (t % disc != 0) return false;
        }
        return true;
    }

    std::deque<std::vector<int>> potential_range_with_length(const int& num,const int& denominator){
        std::deque<std::vector<int>> result;
        std::vector<int> arr;
        int position = num; 
        for (int i = 0; i < num; i++) arr.push_back(0);
        ready_flag:
        result.push_back(arr);
        not_ready_flag:
        arr[position-1] = (arr[position-1] + 1) % denominator;
        if (arr[position-1] != 0){
            position = num;
            goto ready_flag;
        }else{
            position--;
            if (position != 0){
                goto not_ready_flag;
            }else{
                return result;
            }
        }
    }

    bool in_coset_with(const GROUP_ELEMENT& a,const GROUP_ELEMENT& b, const std::deque<GROUP_ELEMENT>& subgroup_generators){
        for (auto coeffs : potential_range_with_length(subgroup_generators.size(),disc)){
            Eigen::Matrix<int,1,N> test_coordinates = a.coordinates - b.coordinates;
            for (int j = 0; j < coeffs.size(); j++) test_coordinates += coeffs[j] * subgroup_generators[j].coordinates;
            bool t = true;
            for (auto cood : test_coordinates){
                if (cood % disc != 0){
                    t = false;
                    break;
                }
            }
            if (t == true) return true;
        }
        return false;
    }

    void generate_dual_group(){
        for (auto item : potential_range_with_length(N,disc)){
            if (integral_check(item) == true){
                dual_group.push_back(GROUP_ELEMENT(item,the_order_of(item)));
            }
        }
    }

    void extract_dual_group_generators(){
        std::deque<GROUP_ELEMENT> remain_list = dual_group;
        dual_group_generators.push_back(*std::max_element(remain_list.begin(),remain_list.end()));
        int products_of_orders = dual_group_generators.front().order;
        while (products_of_orders < disc){
            std::deque<GROUP_ELEMENT> result_list;
            std::deque<std::deque<GROUP_ELEMENT>> cosets;
            for (auto item : remain_list){
                std::deque<GROUP_ELEMENT> new_coset;
                bool coset_found = false;
                for (auto& coset : cosets){
                    if (in_coset_with(item,coset[0],dual_group_generators) == true){
                        coset.push_back(item);
                        coset_found = true;
                        break;
                    }
                }
                if (coset_found == false){
                    new_coset.push_back(item);
                    cosets.push_back(new_coset);
                }
            }
            for (auto coset : cosets) result_list.push_back(*std::min_element(coset.begin(),coset.end()));
            remain_list = result_list;
            if (remain_list.size() != 1){
                dual_group_generators.push_back(*std::max_element(remain_list.begin(),remain_list.end()));
                products_of_orders *= dual_group_generators.back().order;
            }else{
                break;
            }
        }
    }

    public:
    INTEGRAL_LATTICE(Eigen::Matrix<int,N,N> input_matrix){
        dim = N;
        intersection_form = input_matrix;
        disc = intersection_form.determinant();
        generate_dual_group();
        extract_dual_group_generators();
    }

    const Eigen::Matrix<int,N,N>& get_lattice(){return intersection_form;}
    const int& get_dim() {return dim;}
    const int& get_size(){return disc;}

    void get_inform(bool more_info = true){
        std::cout << "The integral lattice with intersection form" << std::endl
                  << intersection_form << std::endl
                  << "with dimension " << dim 
                  << " and discriminant " << disc << std::endl;
        
        if (more_info == true){
        std::cout << "The dual group is of size " << dual_group.size() << " with elements:" << std::endl;
            for (auto term : dual_group)
            {
                std::cout << "[" << term.coordinates << "] is of order " << term.order << std::endl;
            }
        std::cout << "The minimum number of the generators of the dual group is " 
                      << dual_group_generators.size() << " with elements:" << std::endl;
            for (auto term : dual_group_generators){
                std::cout << "[" << term.coordinates << "] is of order " << term.order << std::endl;
            }
        }else{
            std::cout << "The minimum number of the generators of the dual group is " << dual_group_generators.size() << std::endl;
        }
    }
};

int main()
{
    auto start_time = std::chrono::steady_clock::now();
    
    int t = 8;
    
    Eigen::Matrix<int,3,3> A;
    // A << 3, t, 1,
        //  t, 2, 0,
        //  1, 0, 1;
    // A << 3, 1, 0,
        //  1, 7, t,
        //  0, t, 4;
    A << 3, 1, 1,
         1, 7, t,
         1, t, 13;
    INTEGRAL_LATTICE<3> B = A;
    
    // Eigen::Matrix<int,4,4> A;
    // A << 2, 0, 0, 0,
        //  0, 2, 0, 0,
        //  0, 0, 2, 0,
        //  0, 0, 0, 2;
    // INTEGRAL_LATTICE<4> B = A;
    
    B.get_inform();
    
    auto end_time = std::chrono::steady_clock::now();
    double duration_sec = std::chrono::duration<double>(end_time - start_time).count();
    std::cout << "Total time cost: " << duration_sec << " secs." << std::endl;
    std::cin.get();
}