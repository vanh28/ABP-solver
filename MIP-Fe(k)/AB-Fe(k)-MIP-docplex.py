import time
import argparse
import math
import psutil
import threading
from memory_profiler import memory_usage
from docplex.mp.model import Model
from cplex.exceptions import CplexSolverError
from func_timeout import func_timeout, FunctionTimedOut

def calculate_label_differences(labels, edges):
  label_differences = []
  for edge in edges:
    v1, v2 = edge
    label_difference = abs(labels[v1-1] - labels[v2-1])  
    label_differences.append(label_difference)
  return label_differences

class ContextTimer(object):
    def __init__(self, msg):
        self.msg = msg
        self.start = 0
        self.process = psutil.Process()
        self.start_memory = 0
        self.running = True
        
    def __enter__(self):
        self.start = time.time()
        self.start_memory = self.process.memory_info().rss
        print('--> begin {0}'.format(self.msg))
        self.thread = threading.Thread(target=self.print_time_elapsed)
        self.thread.start()
        return self  # return value is value of with ()
        
    def __exit__(self, *args):
        self.running = False
        self.thread.join()
        elapsed = time.time() - self.start
        self.msecs = math.ceil(1000 * elapsed)
        print('<-- end {0},  time: {1:.0f} ms'.format(self.msg, self.msecs))

    def print_time_elapsed(self):
        while self.running:
            elapsed = time.time() - self.start
            current_memory = self.process.memory_info().rss
            memory_used = current_memory - self.start_memory
            memory_used_mb = memory_used / (1024 * 1024)
            if memory_used_mb < 0:
                memory_used_mb = 0  # Đặt giá trị bộ nhớ âm thành 0
            print(f"Thời gian đã trôi qua: {elapsed:.2f} giây, Bộ nhớ đã sử dụng: {memory_used_mb:.2f} MB", end='\r')
            time.sleep(1)

def build_model(vertices, edges, k):
    """
    Xây dựng mô hình MIP cho bài toán antibandwidth theo dạng feasibility model (FE(k))

    Args:
        vertices: Danh sách các đỉnh
        edges: Danh sách các cạnh (cặp đỉnh)
        k: Khoảng cách tối thiểu giữa các nhãn

    Returns:
        Mô hình MIP
    """
    n = len(vertices)  # Số đỉnh
    model = Model(name="Antibandwidth_Feasibility", log_output=True)
    model.parameters.timelimit = 1800 
    model.parameters.workmem = 122880
    print(model.parameters.print_info_to_string())
    # Biến quyết định
    x = model.binary_var_matrix(n, n, name='x')

    # Ràng buộc 1: Mỗi đỉnh được gán đúng một nhãn
    for i in vertices:
        model.add_constraint(model.sum(x[i-1, l-1] for l in vertices) == 1, f"VERTICES_{i}")

    # Ràng buộc 2: Mỗi nhãn được gán cho đúng một đỉnh
    for l in range(n):
        model.add_constraint(model.sum(x[i-1, l] for i in vertices) == 1, f"LABELS_{l}")

    # Ràng buộc OBJ-k: Không được gán nhãn quá gần nhau cho các đỉnh kề nhau

    # Tính độ đo đỉnh
    degree = {i: 0 for i in vertices}
    for (i, j) in edges:
        degree[i] += 1
        degree[j] += 1

    # Tìm đỉnh có độ đo lớn nhất
    max_degree_vertex = max(degree, key=degree.get)

    # Thêm ràng buộc phá vỡ đối xứng
    for l in range((n // 2) + 1, n + 1):
        # Xác định biến tương ứng
        indices = [x[max_degree_vertex-1, l-1]]
        coefficients = [1]
        model.add_constraint(
            model.scal_prod(indices, coefficients) == 0,
            f"SYMMETRY_BREAK_{max_degree_vertex}_{l}"
        )
        
    # Ràng buộc OBJ-k: Không được gán nhãn quá gần nhau cho các đỉnh kề nhau
    model.add_constraints(
        (model.scal_prod(
            [x[i-1, l] for l in range(l2, l2 + k + 1)] + [x[j-1, l] for l in range(l2, l2 + k + 1)], 
            [1] * (2 * (k + 1))
        ) <= 1 for (i, j) in edges for l2 in range(n - k)),
        (f"NO_CLOSE_LABELS_{i}_{j}_{l2}" for (i, j) in edges for l2 in range(n - k))
    )
            
    return model

def solve_antibandwidth(vertices, edges, UB, LB, time_limit=1800):
    """
    Giải bài toán antibandwidth sử dụng mô hình (FE(k))

    Args:
        vertices: Danh sách các đỉnh
        edges: Danh sách các cạnh (cặp đỉnh)

    Returns:
        Giá trị k lớn nhất thỏa mãn (FE(k))
    """
    k = 1
    
    # TÌm giá trị trong khoảng [LB, UB]
    # k = LB - 1
    
    labels = [None] * len(vertices)
    while True:
        # TÌm giá trị trong khoảng [LB, UB]
        # if k >= UB:
        #     break
        
        model = build_model(vertices, edges, k)
        model.parameters.timelimit.set(time_limit)
        try:
            solution = model.solve()
        except CplexSolverError as e:
            print(f"Lỗi khi giải quyết mô hình: {e}")
            break
        
        if solution is None:
            break
        
        for i in range(len(vertices)):
            for l in range(len(vertices)):
                if solution.get_value(model.get_var_by_name(f'x_{i}_{l}')) > 0.5:
                    labels[i] = l + 1  
                    break
        
        print(f"Giá trị Antibandwidth: {k+1}")

        print("Các nhãn:", labels)
        
        # Kiểm tra xem có phải bài toán antibandwidth
        label_differences = calculate_label_differences(labels, edges)
        if any(element < k+1 for element in label_differences):
            print("Không phải Antibandwidth")
        else :
            print("Bài toán Antibandwidth")
        
        k += 1
    return k , labels

def read_mtx_rnd_file(file_path):
    edges = set()

    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Bỏ qua dòng đầu tiên mô tả
    size_info = lines[1].strip().split()
    num_vertices = int(size_info[0])  # Số đỉnh

    # Đọc các cạnh từ các dòng tiếp theo
    for line in lines[2:]:
        row, col = map(int, line.strip().split())
        edges.add((row, col))
    
    return num_vertices, list(edges)

def run_with_timeout_and_memory(vertices, edges, timeout, UB, LB):
    try:
         with ContextTimer("Time and memory usage"):
            # Đo bộ nhớ sử dụng của hàm solve_antibandwidth với giới hạn thời gian
            mem_usage = memory_usage((func_timeout, (timeout, solve_antibandwidth), {'args': (vertices, edges, UB, LB)}))
            print(f"Bộ nhớ sử dụng: {max(mem_usage) - min(mem_usage):.2f} MB")
    except FunctionTimedOut:
        print(f"Hàm solve_antibandwidth đã vượt quá thời gian {timeout} giây")

def main():
    
    #Đọc dữ liệu đầu vào từ tệp .mtx.rnd
    parser = argparse.ArgumentParser(description='Giải bài toán antibandwidth trên đồ thị.')
    parser.add_argument('input_file', type=str, help='Đường dẫn đến tệp đầu vào chứa dữ liệu đồ thị (.mtx.rnd).')
    args = parser.parse_args()
    
    num_vertices, edges, UB, LB = read_mtx_rnd_file(args.input_file)
    vertices = list(range(1, num_vertices+1))  # Danh sách các đỉnh từ 1 đến 48
    timeout = 1800  # Thời gian giới hạn 1800 giây
    
    print(UB)
    print(LB)
    
    run_with_timeout_and_memory(vertices, edges, timeout, UB, LB)
    
    

if __name__ == "__main__":
    main()
