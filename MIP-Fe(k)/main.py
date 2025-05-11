import cplex
from cplex.exceptions import CplexSolverError
import time
import threading
import argparse
import psutil
from memory_profiler import memory_usage
import math
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

def build_model(vertices, edges, k, time_limit=1800):
    """
    Xây dựng mô hình MIP cho bài toán antibandwidth theo dạng (FE(k))

    Args:
        vertices: Danh sách các đỉnh
        edges: Danh sách các cạnh (cặp đỉnh)
        k: Khoảng cách tối thiểu giữa các nhãn

    Returns:
        Mô hình MIP dưới dạng đối tượng cplex.Cplex
    """
    # Số đỉnh
    n = len(vertices) 
    
    # Tạo tham số cho mô hình
    model = cplex.Cplex()
    model.parameters.timelimit.set(time_limit)
    model.parameters.clocktype.set(2)
    model.parameters.workmem.set(122880)
    model.parameters.simplex.pgradient.set(4)
    model.parameters.mip.strategy.file.set(3)
    
    # Tạo biến quyết định
    var_names = [f"x_{i}_{l}" for i in vertices for l in vertices]
    model.variables.add(names=var_names, types=['B'] * len(var_names))

    rows = []
    rhs = []
    senses = []
    names = []
    
    # Ràng buộc 1: Mỗi đỉnh được gán đúng một nhãn
    for i in vertices:
        indices = [var_names[(i-1)*n + l-1] for l in vertices]
        coefficients = [1] * n
        rows.append((indices, coefficients))
        rhs.append(1)
        senses.append('E')
        names.append(f"VERTICES_{i}")
    
    # Ràng buộc 2: Mỗi nhãn được gán cho đúng một đỉnh
    for l in range(1, n + 1):
        indices = [var_names[(i-1)*n + l-1] for i in vertices]
        coefficients = [1] * n
        rows.append((indices, coefficients))
        rhs.append(1)
        senses.append('E')
        names.append(f"LABELS_{l}")
    
    # Ràng buộc OBJ-k: Không được gán nhãn quá gần nhau cho các đỉnh kề nhau
    for (i, j) in edges:
        for l2 in range(n - k):
            indices = [var_names[(i-1)*n + l] for l in range(l2, l2 + k + 1)] + \
                    [var_names[(j-1)*n + l] for l in range(l2, l2 + k + 1)]
            coefficients = [1] * (2 * (k + 1))
            rows.append((indices, coefficients))
            rhs.append(1)
            senses.append('L')
            names.append(f"NO_CLOSE_LABELS_{i}_{j}_{l2}")
            
    # Tính độ đo của các đỉnh
    degree = {i: 0 for i in vertices}
    for (i, j) in edges:
        degree[i] += 1
        degree[j] += 1

    # Tìm đỉnh có độ đo lớn nhất
    max_degree_vertex = max(degree, key=degree.get)
    
    # Ràng buộc phá vỡ đối xứng
    for l in range((n // 2) + 1, n + 1):
        indices = [var_names[(max_degree_vertex-1)*n + l-1]]
        coefficients = [1]
        rows.append((indices, coefficients))
        rhs.append(0)
        senses.append('L')
        names.append(f"SYMMETRY_BREAK_{max_degree_vertex}_{l}")

    # Add constraints in batch
    model.linear_constraints.add(
        lin_expr=rows,
        senses=senses,
        rhs=rhs,
        names=names
    )

    return model



def solve_antibandwidth(vertices, edges, UB, LB):
    """
    Giải bài toán antibandwidth sử dụng mô hình (FE(k))

    Args:
        vertices: Danh sách các đỉnh
        edges: Danh sách các cạnh (cặp đỉnh)

    Returns:
        Giá trị k lớn nhất thỏa mãn (FE(k))
    """
    k = LB - 1
    # k = 1
    labels = [None] * len(vertices)
    while True:
        if k >= UB:
            break
        
        model = build_model(vertices, edges, k)
        
        try:
            model.solve()
            solution = model.solution
            status = solution.get_status_string()
        except CplexSolverError as e:
            print(f"Lỗi: {e}")
            break
        
        print("Trạng thái giải pháp:", status)
        
        if not solution.is_primal_feasible():
            break
        
        var_values = solution.get_values()
        for i in range(len(vertices)):
            for l in range(len(vertices)):
                if var_values[i * len(vertices) + l] > 0.5:
                    labels[i] = l + 1  
                    break
                
        print(f"Giá trị Antibandwidth: {k+1}")
        print(f"Giá trị Antibandwidth: {k+1}")
        print("Các nhãn:", labels)
        
        # Kiểm tra xem có phải bài toán antibandwidth
        label_differences = calculate_label_differences(labels, edges)
        if any(element < k + 1 for element in label_differences):
            print("Không phải bài toán Antibandwidth")
        else :
            print("Bài toán Antibandwidth")
        k += 1
    return k, labels

def read_mtx_rnd_file(file_path):
    edges = set()

    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Bỏ qua dòng đầu tiên mô tả
    size_info = lines[1].strip().split()
    num_vertices = int(size_info[0])  # Số đỉnh
    LB = int(size_info[3])  # Upper bound
    UB = int(size_info[4])  # Lower bound
    # Đọc các cạnh từ các dòng tiếp theo
    for line in lines[2:]:
        row, col = map(int, line.strip().split())
        edges.add((row, col))
    
    return num_vertices, list(edges), UB, LB

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
