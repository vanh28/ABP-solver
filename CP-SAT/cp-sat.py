from ortools.sat.python import cp_model
import argparse

def calculate_label_differences(labels, edges):
  label_differences = []
  for edge in edges:
    v1, v2 = edge
    # print(v1, v2, labels[v1-1], labels[v2-1], abs(labels[v1-1] - labels[v2-1]))
    label_difference = abs(labels[v1-1] - labels[v2-1])  
    label_differences.append(label_difference)
  return label_differences

def solveABP(num_vertices, edges, LB, UB):
     # Tạo model
    model = cp_model.CpModel()
    
    # Biến quyết định
    labels = [model.NewIntVar(0, num_vertices - 1, f'label_{i}') for i in range(0, num_vertices)]
    b = model.NewIntVar(LB, UB, 'b')
    
    # Ràng buộc C.2: b ≤ abs(labels[u] - labels[v])
    for u, v in edges:
        diff = model.NewIntVar(1, num_vertices, f'diff_{u}_{v}')
        model.AddAbsEquality(diff, labels[u - 1] - labels[v - 1])
        model.Add(b <= diff)
    
    # Ràng buộc C.3: tất cả nhãn phải khác nhau
    model.AddAllDifferent(labels)
    
    # Tính bậc của từng đỉnh và tìm đỉnh có bậc cao nhất
    vertex_degree = [0] * num_vertices
    for u, v in edges:
        vertex_degree[u - 1] += 1
        vertex_degree[v - 1] += 1
    
    max_degree_vertex = max(range(num_vertices), key=lambda v: vertex_degree[v])
    
    # Ràng buộc phá vỡ đối xứng
    model.Add(labels[max_degree_vertex] <= num_vertices // 2)
    
    # Hàm mục tiêu: tối đa hóa b
    model.Maximize(b)
    
    # Giải bài toán
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    solver.parameters.max_time_in_seconds = 2
    # Xuất kết quả
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        print(f"Objective value (b): {solver.Value(b)}")
        print("Labels:", [solver.Value(label) for label in labels])
        calculate_diff = calculate_label_differences([solver.Value(label) for label in labels], edges)
        # neu calculate_diff co gia tri bang > b thi khong phai AB
        print("Label differences:", calculate_diff)
        
    else:
        print("No feasible solution found.")

def read_mtx_rnd_file(file_path):
    edges = set()

    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    size_info = lines[1].strip().split()
    num_vertices = int(size_info[0]) 
    LB = int(size_info[3]) 
    UB = int(size_info[4])
    
    for line in lines[2:]:
        row, col = map(int, line.strip().split())
        edges.add((row, col))
    
    return num_vertices, list(edges), LB, UB

def main():
    parser = argparse.ArgumentParser(description='Giải bài toán antibandwidth')
    parser.add_argument('input_file', type=str, help='Đường dẫn đến file đầu vào')
    args = parser.parse_args()

    
    num_vertices, edges, LB, UB = read_mtx_rnd_file(args.input_file)
    
    
    print(f"Number of vertices: {num_vertices}")
    print(f"Number of edges: {len(edges)}")
    print(f"Upper Bound (UB): {UB}")
    print(f"Lower Bound (LB): {LB}")
    
    solveABP(num_vertices, edges, LB, UB)
    
if __name__ == "__main__":
    main()
