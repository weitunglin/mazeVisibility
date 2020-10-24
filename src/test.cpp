#include <iostream>
using namespace std;

bool get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y, 
    double p2_x, double p2_y, double p3_x, double p3_y, double *i_x, double *i_y) {
    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0)
    {
        // Collision detected
        if (i_x != 0)
            *i_x = p0_x + (t * s1_x);
        if (i_y != 0)
            *i_y = p0_y + (t * s1_y);
        return true;
    }

    return false; // No collision
}

int main() {
    double s[2];
    bool result = get_line_intersection(6.80725, 2.6609, 31.0075, 20.3907, 10.0, 5.0, 15.0, 5.0, &s[0], &s[1]);

    cout << result << endl;
    cout << s[0] << " " << s[1] << endl;
}
