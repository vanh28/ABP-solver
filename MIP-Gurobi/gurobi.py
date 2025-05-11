from gurobipy import Model, GRB
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
        label_difference = abs(labels[v1 - 1] - labels[v2 - 1])
        label_differences.append(label_difference)
    return label_differences

def build_model(vertices, edges, k, time_limit=1800):
    
    n = len(vertices)
    model = Model("Antibandwidth")
    model.setParam('TimeLimit', time_limit)

    x = model.addVars(vertices, vertices, vtype=GRB.BINARY, name="x")

    # Ràng buộc 1: Mỗi đỉnh được gán đúng một nhãn
    for i in vertices:
        model.addConstr(sum(x[i, l] for l in vertices) == 1, name=f"VERTICES_{i}")

    # Ràng buộc 2: Mỗi nhãn được gán đúng một đỉnh
    for l in vertices:
        model.addConstr(sum(x[i, l] for i in vertices) == 1, name=f"LABELS_{l}")

    # Ràng buộc OBJ-k: Không được gán nhãn quá gần nhau cho các đỉnh kề nhau
    for (i, j) in edges:
        for l2 in range(1, n - k + 1):
            model.addConstr(
                sum(x[i, l] for l in range(l2, l2 + k + 1)) + 
                sum(x[j, l] for l in range(l2, l2 + k + 1)) <= 1, 
                name=f"NO_CLOSE_LABELS_{i}_{j}_{l2}"
            )

    # Ràng buộc phá vỡ đối xứng
    degree = {i: 0 for i in vertices}
    for (i, j) in edges:
        degree[i] += 1
        degree[j] += 1
    max_degree_vertex = max(degree, key=degree.get)
    for l in range((n // 2) + 1, n + 1):
        model.addConstr(x[max_degree_vertex, l] == 0, name=f"SYMMETRY_BREAK_{max_degree_vertex}_{l}")

    return model

def solve_antibandwidth(vertices, edges, UB, LB):
    k = LB - 1
    labels = [None] * len(vertices)
    while True:
        if k >= UB:
            break
        
        model = build_model(vertices, edges, k)
        model.optimize()

        if model.Status != GRB.OPTIMAL:
            print(f"Không tìm thấy lời giải với giá trị AB = {k + 1}.")
            break

        var_values = model.getAttr('X', model.getVars())
        for i in range(len(vertices)):
            for l in range(len(vertices)):
                if var_values[i * len(vertices) + l] > 0.5:
                    labels[i] = l + 1
                    break

        print(f"Giá trị Antibandwidth: {k+1}")
        print("Các nhãn:", labels)

        label_differences = calculate_label_differences(labels, edges)
        if any(element < k + 1 for element in label_differences):
            print("Không phải bài toán Antibandwidth")
        else:
            print("Bài toán Antibandwidth")
        k += 1
    return k, labels

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
    
    return num_vertices, list(edges), UB, LB

def main():
    parser = argparse.ArgumentParser(description='Giải bài toán antibandwidth trên đồ thị.')
    parser.add_argument('input_file', type=str, help='Đường dẫn đến tệp đầu vào chứa dữ liệu đồ thị (.mtx.rnd).')
    args = parser.parse_args()

    num_vertices, edges, UB, LB = read_mtx_rnd_file(args.input_file)
    vertices = list(range(1, num_vertices + 1))

    solve_antibandwidth(vertices, edges, UB, LB)

if __name__ == "__main__":
    main()
