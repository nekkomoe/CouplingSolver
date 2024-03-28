/**
 * 解决非结构网格的计算问题
*/

// 节点的物理场值
#define NodeFieldContent \
    Ts, Tc, Phi, P
// 边的物理场值
#define EdgeFieldContent \
    U

enum class NodeField { NodeFieldContent, _size};
enum class EdgeField { EdgeFieldContent, _size};
constexpr int NodeFieldNum = static_cast<int>(NodeField::_size);
constexpr int EdgeFieldNum = static_cast<int>(EdgeField::_size);

union NodeFieldValue {
    // 体积元重心处的值
    double value[NodeFieldNum];
    struct { double NodeFieldContent; };
};
union EdgeFieldValue {
    // 交错网格（边上值）
    double value[EdgeFieldNum];
    struct { double EdgeFieldContent; };
};

// 声明后面会用到的结构体
struct MeshNode;
struct MeshEdge;


struct MeshEdge {
    bool is_boundary;              // 是否是边界
    MeshNode *from_node, *to_node; // 边的两个节点
    double from_dis, to_dis;       // 体积元到边的垂直距离
    double x1, y1, x2, y2;         // 边的坐标
    double dL;                     // 边的长度
    EdgeFieldValue _edge_value;    // 边的物理场值
    NodeFieldValue _node_value;    // 体积元在边上的值
    NodeFieldValue _node_grad;     // 体积元在边上的梯度
    EdgeFieldValue F_value;        // 非线性函数值
    void clear() {
        from_node = to_node = nullptr;
        from_dis = to_dis = 0;
        is_boundary = false;
    }
    void geoinfo();                // 计算边的几何信息
    void interp();                 // 计算体积元在边上的值
    void grad();                   // 计算体积元的梯度
    void grad_upwind(int U_idx);   // 计算体积元的梯度（迎风格式）
};

struct MeshNode {
    typedef pair<MeshEdge*, MeshNode*> edge_t;
    bool is_boundary;              // 是否是边界
    int node_idx;                  // 节点编号
    double x, y;                   // 节点坐标
    double dV;                     // 体积元的体积
    NodeFieldValue _node_value;    // 节点的物理场值
    NodeFieldValue F_value;        // 非线性函数值
    vector<edge_t> neighbor;       // 邻居节点
    void gemoinfo() {
        // 计算体积元的体积
        dV = 0;
        for(auto &e : neighbor) {
            auto [edge, to_node] = e;
            double x1 = edge->x1-x, y1 = edge->y1-y,
                   x2 = edge->x2-x, y2 = edge->y2-y;
            // 计算三角面积, 假设体积元都是凸多边形
            double tri = fabs(x1*y2 - x2*y1);
            dV += tri;
        }
    }
};

struct Mesh {
    vector<MeshNode> node, bc_node;
    vector<MeshEdge> edge, bc_edge;
    void init();   // 初始化网格
    void update(); // 更新网格信息
};


double point2line(double x, double y, double x1, double y1, double x2, double y2) {
    // 计算点到直线的距离
    double A = y2 - y1;
    double B = x1 - x2;
    double C = x2 * y1 - x1 * y2;
    return fabs(A*x + B*y + C) / sqrt(A*A + B*B);
}

// MeshEdge 函数声明
void MeshEdge::geoinfo() {
    // 计算边的几何信息, 其中(x1,y1)-(x2,y2)为边的坐标
    if(from_node) {
        from_dis = point2line(from_node->x, from_node->y, x1, y1, x2, y2);
    }
    if(to_node) {
        to_dis   = point2line(to_node->x, to_node->y,     x1, y1, x2, y2);
    }
    // 计算是否为边界: 如果from_node或to_node为空, 则是边界
    is_boundary  = (!from_node) || (!to_node);
    // 计算长度
    dL = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}
void MeshEdge::interp() {
    // 计算体积元在边上的值
    for(int i = 0 ; i < NodeFieldNum; i++) {
        double from_value = from_node->_node_value.value[i],
                to_value   = to_node->_node_value.value[i];
        double dis_p      = from_dis/(from_dis+to_dis);
        // 距离越远，权重越小
        double value      = from_value*(1-dis_p) + to_value*dis_p;
        _node_value.value[i] = value;
    }
}
void MeshEdge::grad() {
    // 计算体积元的梯度
    for(int i = 0 ; i < NodeFieldNum; i++) {
        double from_value = from_node->_node_value.value[i],
                to_value   = to_node->_node_value.value[i];
        // 方向向外
        double value      = (to_value - from_value) / (to_dis + from_dis);
        _node_grad.value[i] = value;
    }
}
void MeshEdge::grad_upwind(int U_idx = 0) {
    // 前置要求：interp()已经被调用, 体积元在边上的值已经计算
    // U_idx: 速度场的索引
    // 计算体积元梯度（迎风格式）
    for(int i = 0 ; i < NodeFieldNum; i++) {
        double from_value = from_node->_node_value.value[i],
                to_value   = to_node->_node_value.value[i];
        double U          = _edge_value.value[U_idx];
        double value      = 0;
        if(U > 0) { // 速度向外
            value = (_node_value.value[i] - from_value) / from_dis;
        } else {
            value = (to_value - _node_value.value[i]) / to_dis;
        }
        _node_grad.value[i] = value;
    }
}

// Mesh 函数声明
void Mesh::init() {
    // 初始化边的连接信息
    for(auto &e : edge)    { e.clear(); }
    for(auto &e : bc_edge) { e.clear(); }
    for(auto &n : node) {
        for(auto &e : n.neighbor) {
            // 这里会被计算两次, 最后一次计算时确定方向
            auto [edge, to_node] = e;
            edge->from_node = &n;
            edge->to_node   = to_node;
        }
    }
    // 初始化边的几何信息
    for(auto &e : edge)    { e.geoinfo(); }
    for(auto &e : bc_edge) { e.geoinfo(); }
    // 初始化节点的几何信息
    for(auto &n : node)    { n.gemoinfo(); }
    for(auto &n : bc_node) { n.gemoinfo(); }
}
void Mesh::update() {
    // 更新物理场信息
    // 1) 计算边上的插值, 梯度, 迎风梯度
    for(auto &e: edge) { e.interp(); }
    for(auto &e: edge) { e.grad(); e.grad_upwind(); }
}

// 生成Jacobian的非零元素结构
void umJacNNZ() {
}

// 计算非结构网格的非线性函数
void umF() {

}

void unit_test() {
}
