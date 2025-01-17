#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <limits>
#include <map>
#include <queue>
#include <fstream>

// Константы
const int NUM_POINTS = 100; // Количество точек
const double RADIUS = 100.0; // Радиус окружности
const double DISTANCE_COST = 10.0; // Стоимость единицы расстояния
// Число ПИ
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

struct Point {
    double x, y; // Координаты точки
    int id; // Уникальный ID точки
};

// Генерация случайных точек внутри окружности
std::vector<Point> generatePoints(int numPoints, double radius) {
    std::vector<Point> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    std::uniform_real_distribution<> angleDist(0, 2 * M_PI);

    for (int i = 0; i < numPoints; ++i) {
        double r = radius * std::sqrt(dist(gen));
        double theta = angleDist(gen);
        points.push_back({r * std::cos(theta), r * std::sin(theta), i});
    }
    return points;
}

// Вычисление расстояния между двумя точками
double distance(const Point& a, const Point& b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

// Построение графа: соединяем каждую точку с 2-6 ближайшими соседями
std::map<int, std::vector<std::pair<int, double>>>
buildGraph(const std::vector<Point>& points) {
    std::map<int, std::vector<std::pair<int, double>>> graph;

    for (const auto& point : points) {
        std::vector<std::pair<double, int>> distances;

        for (const auto& other : points) {
            if (point.id != other.id) {
                distances.push_back({distance(point, other), other.id});
            }
        }

        // Сортируем по расстоянию
        std::sort(distances.begin(), distances.end());

        // Соединяем с 2-6 ближайшими соседями
        int numConnections = 2 + rand() % 5; // От 2 до 6
        for (int i = 0; i < numConnections && i < distances.size(); ++i) {
            graph[point.id].push_back({distances[i].second, distances[i].first});
        }
    }

    return graph;
}

// Функция для поиска кратчайшего пути с использованием алгоритма Дейкстры
std::vector<int> findShortestPath(const std::map<int, std::vector<std::pair<int, double>>>& graph, int start, int end) {
    std::map<int, double> minDistance;
    std::map<int, int> previous;
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;

    // Инициализация
    for (const auto& [node, _] : graph) {
        minDistance[node] = std::numeric_limits<double>::max();
    }
    minDistance[start] = 0;
    pq.push({0, start});

    while (!pq.empty()) {
        auto [currentDistance, currentNode] = pq.top();
        pq.pop();

        // Если достигли конечной точки
        if (currentNode == end) {
            break;
        }

        for (const auto& [neighbor, cost] : graph.at(currentNode)) {
            double newDistance = currentDistance + cost;
            if (newDistance < minDistance[neighbor]) {
                minDistance[neighbor] = newDistance;
                previous[neighbor] = currentNode;
                pq.push({newDistance, neighbor});
            }
        }
    }

    // Восстановление пути
    std::vector<int> path;
    for (int at = end; at != start; at = previous[at]) {
        path.push_back(at);
    }
    path.push_back(start);
    std::reverse(path.begin(), path.end());

    return path;
}

// Вывод графа в файл для визуализации
void visualizeGraph(const std::vector<Point>& points,
                    const std::map<int, std::vector<std::pair<int, double>>>& graph,
                    const std::vector<int>& path) {
    std::ofstream file("graph.dot");
    file << "graph G {\n";

    // Рисуем рёбра
    for (const auto& [node, neighbors] : graph) {
        for (const auto& [neighbor, cost] : neighbors) {
            file << "  " << node << " -- " << neighbor << " [label=\"" << cost << "\"];\n";
        }
    }

    // Подсвечиваем путь
    for (size_t i = 0; i < path.size() - 1; ++i) {
        file << "  " << path[i] << " -- " << path[i + 1] << " [color=red, penwidth=2.0];\n";
    }

    file << "}\n";
    file.close();
}

int main() {
    // Генерация случайных точек
    auto points = generatePoints(NUM_POINTS, RADIUS);

    // Построение графа
    auto graph = buildGraph(points);

    // Ввод стартовой и конечной точки
    int startPoint, endPoint;
    std::cout << "Enter the start point (0 to " << NUM_POINTS - 1 << "): "; 
    std::cin >> startPoint; 
    std::cout << "Enter the end point (0 to " << NUM_POINTS - 1 << "): "; 
    std::cin >> endPoint;

    // Поиск кратчайшего пути
    auto path = findShortestPath(graph, startPoint, endPoint);

    // Визуализация графа
    visualizeGraph(points, graph, path);

    std::cout << "Path: ";
    for (auto node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    return 0;
}