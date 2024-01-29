// require latest C++ compling standard
#include <iostream>
#include <iomanip>
#include <string>
#include <array>
#include <vector>
#include <deque>
#include <cmath>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <initializer_list>

template<typename T>
class INTEGRAL_LATTICE{
    protected:
    class group_element{
        public:
        std::vector<T> element;
        int order;
        const bool operator<(const group_element& other) const{return this->order < other.order;}
        const bool operator>(const group_element& other) const{return this->order > other.order;}
    };
    private:
    std::vector<std::vector<T>> intersection_matrix;
    int rank;
    int det;
    std::vector<group_element> group_list;
    std::vector<group_element> generator_list;
    const std::vector<std::vector<T>> sub_matrix(const int& i, const std::vector<std::vector<T>>& input_matrix) const{
        std::vector<std::vector<T>> result_matrix;
        result_matrix.reserve(input_matrix.size()-1);
        for (int m = 0; m < input_matrix.size()-1; m++){
            std::vector<int> temp;
            temp.reserve(input_matrix.size()-1);
            for (int n = 0; n < input_matrix.size()-1; n++){
                if (n < i) temp.emplace_back(input_matrix[m+1][n]);
                else       temp.emplace_back(input_matrix[m+1][n+1]);
            }
            result_matrix.emplace_back(temp);
        }
        return result_matrix;
    }
    const int determinant(const std::vector<std::vector<T>>& input_matrix) const{
        int rank = input_matrix.size();
        int result = 0;
        if (rank == 1) return input_matrix[0][0];
        else if (rank == 2) return input_matrix[0][0]*input_matrix[1][1] - input_matrix[0][1]*input_matrix[1][0];
        else if (rank == 3){
            result += input_matrix[0][0] * input_matrix[1][1] * input_matrix[2][2];
            result -= input_matrix[0][0] * input_matrix[2][1] * input_matrix[1][2];
            result -= input_matrix[1][0] * input_matrix[0][1] * input_matrix[2][2];
            result += input_matrix[1][0] * input_matrix[2][1] * input_matrix[0][2];
            result += input_matrix[2][0] * input_matrix[0][1] * input_matrix[1][2];
            result -= input_matrix[2][0] * input_matrix[1][1] * input_matrix[0][2];
            return result;
        }else{
            for (int i = 0; i < input_matrix.size(); i++){
                if (i % 2 == 0) result += input_matrix[0][i] * determinant(sub_matrix(i,input_matrix));
                else            result -= input_matrix[0][i] * determinant(sub_matrix(i,input_matrix));
            }
            return result;
        }
    }
    const int the_order_of(const std::vector<T>& input_vector) const{
        for (int i = 1; i <= abs(det); i++){
            int t = 0;
            for (int k = 0; k < rank; k++){
                if ((input_vector[k] * i) % det != 0) break;
                else t++;
            }
            if (t == rank) return i;
        }
        return 0;
    }
    template<typename M>
    const bool in_integral_dual(const M& input_vector) const{
        for (int i = 0; i < rank; i++){
            int t = 0;
            for (int j = 0; j < rank; j++) t += input_vector[j] * intersection_matrix[i][j];
            if (t % det != 0) return false;
        }
        return true;
    }
    std::vector<group_element> group_list_checker(){
        int length = rank;
        int lower_bound = 0;
        int upper_bound = abs(det);
        std::vector<T> steps(length,1);
        for (int i = 0; i < length; i++){
            steps[i] = intersection_matrix[i][0];
            for (int j = 1; j < length; j++){
                steps[i] = std::gcd(steps[i], intersection_matrix[i][j]);
            }
        }
        std::vector<group_element> result;
        result.reserve(abs(det));
        std::vector<T> arr(length,lower_bound);
        int position = length-1;
        ready_flag:
        if (in_integral_dual(arr)) result.emplace_back(group_element(arr,the_order_of(arr)));
        if (result.size() == abs(det)) return result;
        not_ready_flag:
        arr[position] += steps[position];
        if (arr[position] < upper_bound){
            position = length-1;
            goto ready_flag;
        }else{
            arr[position] = lower_bound;
            position--;
            if (position != -1) goto not_ready_flag;
            else return result;
        }
    }
    const bool in_coset_with(const group_element& a, const group_element& b,
                            const std::vector<group_element>& subgroup_generators) const{
        int length = subgroup_generators.size();
        int lower_bound = 0;
        int upper_bound = abs(det);
        std::vector<T> arr(length,lower_bound);
        int position = length-1;
        bool t;
        ready_flag:
        t = true;
        std::vector<T> test_coordinates;
        test_coordinates.reserve(rank);
        for (int i = 0; i < rank; i++) test_coordinates.push_back(a.element[i] - b.element[i]);
        for (int j = 0; j < length; j++){
            for (int k = 0; k < rank; k++) test_coordinates[k] = test_coordinates[k] + arr[j] * ((subgroup_generators[j].element)[k]);
        }
        for (auto cood : test_coordinates){
            if (cood % det != 0){
                t = false;
                break;
            }
        }
        if (t == true) return true;
        not_ready_flag:
        arr[position] = arr[position]+1;
        if (arr[position] != upper_bound){
            position = length-1;
            goto ready_flag;
        }else{
            arr[position] = lower_bound;
            position--;
            if (position != -1) goto not_ready_flag;
            else return false;
        }
    }
    std::vector<group_element> generators_extractor(){
        std::vector<group_element> found_generators;
        std::vector<group_element> remain_list = group_list;
        found_generators.push_back(*std::max_element(remain_list.begin(),remain_list.end()));
        int products_of_orders = found_generators.front().order;
        while (products_of_orders < std::abs(det)){
            std::vector<group_element> result_list;
            result_list.reserve(abs(det) / products_of_orders);
            std::vector<std::vector<group_element>> cosets;
            cosets.reserve(abs(det) / products_of_orders);
            for (auto item : remain_list){
                std::vector<group_element> new_coset;
                bool coset_found = false;
                for (auto& coset : cosets){
                    if (coset.size() == products_of_orders) continue;
                    if (in_coset_with(item,coset[0],found_generators) == true){
                        coset.emplace_back(item);
                        coset_found = true;
                        break;
                    }
                }
                if (coset_found == false){
                    new_coset.push_back(item);
                    cosets.emplace_back(new_coset);
                    cosets.back().reserve(products_of_orders);
                }
            }
            for (auto coset : cosets) result_list.emplace_back(*std::min_element(coset.begin(),coset.end()));
            remain_list = result_list;
            if (remain_list.size() != 1){
                found_generators.push_back(*std::max_element(remain_list.begin(),remain_list.end()));
                products_of_orders *= found_generators.back().order;
            }else break;
        }
        return found_generators;
    }
    std::vector<T> primary_sub(std::vector<std::vector<T>>& input_matrix) const{
        int rk = input_matrix[0].size();
        std::vector<T> result;
        for (int i = 0; i < rk; i++){
            std::vector<std::vector<T>> primary_sub;
            primary_sub.reserve(i);
            for (int j = 0; j <= i; j++){
                std::vector<T> temp;
                temp.reserve(i);
                for (int k = 0; k <= i; k++) temp.push_back(input_matrix[j][k]);
                primary_sub.push_back(temp);
            }
            result.push_back(determinant(primary_sub));
        }
        return result;
    }
    const bool positive_definite(std::vector<std::vector<T>>& input_matrix) const{
        for (auto term : primary_sub(input_matrix)) if (term <= 0) return false;
        return true;
    }
    std::vector<T> all_sub(std::vector<std::vector<T>>& input_matrix) const{
        int length = input_matrix.size();
        int lower_bound = 0;
        int upper_bound = 2;
        std::vector<T> result;
        result.reserve(int(std::pow(2,length)));
        std::vector<T> arr(length,lower_bound);
        int position = length-1;
        ready_flag:
        std::vector<std::vector<T>> sub_matrix;
        for (int i = 0; i < length; i++){
            if (arr[i] == 0) continue;
            std::vector<T> temp;
            for (int j = 0; j < length; j++) if (arr[j] == 1) temp.push_back(input_matrix[i][j]);
            sub_matrix.push_back(temp);
        }
        result.push_back(determinant(sub_matrix));
        not_ready_flag:
        arr[position] = arr[position]+1;
        if (arr[position] != upper_bound){
            position = length-1;
            goto ready_flag;
        }else{
            arr[position] = lower_bound;
            position--;
            if (position != -1) goto not_ready_flag;
            else return result;
        }
    }
    const bool semi_positive_definite(std::vector<std::vector<T>>& input_matrix) const{
        for (auto term : all_sub(input_matrix)) if (term < 0) return false;
        return true;
    }
    public:
    const bool has_root() const{
        std::vector<std::vector<T>> test_matrix = intersection_matrix;
        if (positive_definite(test_matrix) == false) return false;
        int bound = abs(intersection_matrix[0][0]);
        for (int i = 1; i < rank; i++) bound = std::max(bound,abs(intersection_matrix[i][i]));
        int length = rank;
        int lower_bound = -bound;
        int upper_bound = bound+1;
        std::vector<T> arr(length,lower_bound);
        int position = length-1;
        ready_flag:
        if (quadratic_product(arr) == 2) return true;
        not_ready_flag:
        arr[position] = arr[position]+1;
        if (arr[position] != upper_bound){
            position = length-1;
            goto ready_flag;
        }else{
            arr[position] = lower_bound;
            position--;
            if (position != -1) goto not_ready_flag;
            else return false;
        }
        return false;
    }
    public:
    INTEGRAL_LATTICE(std::initializer_list<T> input_list){
        std::vector<T> lists = input_list;
        rank = int(std::sqrt(input_list.size()));
        for (int j = 0; j < rank; j++){
            std::vector<T> temp;
            temp.reserve(rank);
            for (int k = 0; k < rank; k++) temp.push_back(lists[rank*j + k]);
            intersection_matrix.push_back(temp);
        }
        det = determinant(intersection_matrix);
        group_list = group_list_checker();
        generator_list = generators_extractor();
    }
    INTEGRAL_LATTICE(const std::vector<T>& input_list){
        rank = int(std::sqrt(input_list.size()));
        for (int j = 0; j < rank; j++){
            std::vector<T> temp;
            temp.reserve(rank);
            for (int k = 0; k < rank; k++) temp.push_back(input_list[rank*j + k]);
            intersection_matrix.push_back(temp);
        }
        det = determinant(intersection_matrix);
        group_list = group_list_checker();
        generator_list = generators_extractor();
    }
    const int get_rank()const{return rank;}
    const int get_det() const{return det;}
    void get_matrix(){
        for (auto i : intersection_matrix){
            for (auto j : i) std::cout << j << " ";
            std::cout << std::endl;
        }
    }
    const bool operator>(const T& num) const{
        std::vector<std::vector<T>> test_matrix;
        for (int i = 0; i < rank; i++){
            std::vector<T> temp;
            for (int j = 0; j < rank; j++){
                if (i != j) temp.push_back(intersection_matrix[i][j]);
                else temp.push_back(intersection_matrix[i][j] - num);
            }
            test_matrix.push_back(temp);
        }
        return positive_definite(test_matrix);
    }
    const bool operator>=(const T& num) const{
        std::vector<std::vector<T>> test_matrix;
        for (int i = 0; i < rank; i++){
            std::vector<T> temp;
            for (int j = 0; j < rank; j++){
                if (i != j) temp.push_back(intersection_matrix[i][j]);
                else temp.push_back(intersection_matrix[i][j] - num);
            }
            test_matrix.push_back(temp);
        }
        return semi_positive_definite(test_matrix);
    }
    const bool operator<(const T& num) const{
        std::vector<std::vector<T>> test_matrix;
        for (int i = 0; i < rank; i++){
            std::vector<T> temp;
            for (int j = 0; j < rank; j++){
                if (i != j) temp.push_back(-intersection_matrix[i][j]);
                else temp.push_back(-intersection_matrix[i][j] + num);
            }
            test_matrix.push_back(temp);
        }
        return positive_definite(test_matrix);
    }
    const bool operator<=(const T& num) const{
        std::vector<std::vector<T>> test_matrix;
        for (int i = 0; i < rank; i++){
            std::vector<T> temp;
            for (int j = 0; j < rank; j++){
                if (i != j) temp.push_back(-intersection_matrix[i][j]);
                else temp.push_back(-intersection_matrix[i][j] + num);
            }
            test_matrix.push_back(temp);
        }
        return semi_positive_definite(test_matrix);
    }
    const int quadratic_product(const std::vector<T>& input_element) const{
        int result = 0;
        for (int i = 0; i < rank; i++){
            for (int j = 0; j < rank; j++) result += input_element[i] * intersection_matrix[i][j] * input_element[j];
        }
        return result;
    }
    void show_info(){
        std::cout << "The intersection form with the following intersection matrix" << std::endl;
        get_matrix();
        std::cout << "is of rank " << rank << " and determinant (or discriminant) " << det << std::endl;
        auto info_loop = [](std::vector<group_element>& input_list){
            for (auto i : input_list){
                std::cout << "[ ";
                for (auto j : i.element) std::cout << j << " ";
                std::cout << "] with order " << i.order << std::endl;
            }
        };
        std::cout << "whose integral dual is of the same size " << group_list.size() << " and with elements:" << std::endl;
        info_loop(group_list);
        std::cout << "A set of generators may be of size " << generator_list.size() << " with elements:" << std::endl;
        info_loop(generator_list);
    }
};

int main(){
    auto start_time = std::chrono::steady_clock::now();

    INTEGRAL_LATTICE<int> A = { 2, 0, 0, 0, 0,
                                0, 2, 0, 0, 0,
                                0, 0, 2, 0, 0,
                                0, 0, 0, 2, 0,
                                0, 0, 0, 0, 2};
    A.show_info();
    std::cout << (A > -1) << (A > 1) << (A > 2) << (A >= 2) << (A < 2) << (A <= 2) << (A < 3) << std::endl;

    // INTEGRAL_LATTICE<int> B = {2,  0,  0, -1,  0,  0,
    //                            0,  2, -1,  0,  0,  0,
    //                            0, -1,  2, -1,  0,  0,
    //                           -1,  0, -1,  2, -1,  0,
    //                            0,  0,  0, -1,  2, -1,
    //                            0,  0,  0,  0, -1,  2};
    // B.show_info();
    // std::cout << (B > 0) << std::endl;
    // std::cout << (B < 0) << std::endl;
    // std::cout << B.has_root() << std::endl;

    // INTEGRAL_LATTICE<int> C = {2,  0,  0, -1,  0,  0,  0,
    //                            0,  2, -1,  0,  0,  0,  0,
    //                            0, -1,  2, -1,  0,  0,  0,
    //                           -1,  0, -1,  2, -1,  0,  0,
    //                            0,  0,  0, -1,  2, -1,  0,
    //                            0,  0,  0,  0, -1,  2, -1,
    //                            0,  0,  0,  0,  0, -1,  2};
    // C.show_info();

    // int t = 9;
    // INTEGRAL_LATTICE<int> D = {3, 1, 1,
    //                            1, 7, t,
    //                            1, t, 13};
    // D.show_info();
    // std::cout << D.has_root() << std::endl;
    // std::cout << (D < 19) << " " << (D > -1) << std::endl;

    // for(int i = -9; i < 10; i++){
    //     INTEGRAL_LATTICE<int> E = {3, 1, 1,
    //                                1, 7, i,
    //                                1, i, 13};
    //     if ((E > 0) == true && E.has_root() == false){
    //         std::cout << i << " " << std::endl;
    //     }
    // }
    // INTEGRAL_LATTICE<int> F = { 2, 0, 0, 0, 0, 0,
    //                             0, 2, 0, 0, 0, 0,
    //                             0, 0, 2, 0, 0, 0,
    //                             0, 0, 0, 2, 0, 0,
    //                             0, 0, 0, 0, 2, 0,
    //                             0, 0, 0, 0, 0, 2 };
    // F.show_info();

    auto end_time = std::chrono::steady_clock::now();
    double duration_sec = std::chrono::duration<double>(end_time - start_time).count();
    std::cout << "Total time cost: " << duration_sec << " secs." << std::endl;
    std::cin.get();
}