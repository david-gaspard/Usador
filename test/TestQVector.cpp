/****
 * @date Created on 2025-07-01 at 18:29:17 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the QVector structure.
 * Expected values are produced by the Mathematica notebook: ~/Travail/Brouillons/2022-12-27_Nazarov-Method/math/2025-06-27_Check-Usador.nb
 ***/
#include "QVector.hpp"

/**
 * Compute the relative error (in absolute value) between two qvectors.
 */
double errorQVector(QVector a, QVector a_expc) {
    return abs((a - a_expc).norm())/abs(a.norm());
}

/**
 * Test the basic operators of qvectors (addition, dot product, cross product).
 */
int basicOperations() {
    std::cout << "====== TEST BASIC OPERATIONS ======\n";
    
    QVector a = QVector(1.+2.*I, dcomplex(3., 4.), 5.+6.*I);
    std::cout << TAG_INFO << "a = " << a << "\n";
    
    // Test the conversion of angle parameters (theta, eta) to vector form:
    QVector b = QVector(0.1 + 0.5*I, -0.5 + 0.3*I);
    QVector b_expc = QVector(dcomplex(-0.54839735782824772,0.32591500263646511), 
                             dcomplex(-0.11257475280623035,-0.51849199947310398),
                             dcomplex(1.03687829202776305,0.11608095045981376));
    
    std::cout << TAG_INFO << "b = " << b << ", error = " << errorQVector(b, b_expc) << "\n";
    
    // Test the norm:
    std::cout << TAG_INFO << "Norm b = " << b.norm() << " (expected = 1)\n";
    
    // Test the addition:
    QVector s = a + b;
    QVector s_expc = QVector(dcomplex(0.4516026421717523,2.325915002636465),
                             dcomplex(2.88742524719377,3.481508000526896),
                             dcomplex(6.0368782920277635,6.116080950459814));
    std::cout << TAG_INFO << "Addition:      a + b = " << s << ", error = " << errorQVector(s, s_expc) << "\n";
    
    // Test the subtraction:
    QVector m = a - b;
    QVector m_expc = QVector(dcomplex(1.5483973578282477,1.674084997363535),
                             dcomplex(3.11257475280623,4.518491999473104),
                             dcomplex(3.963121707972237,5.883919049540186));
    std::cout << TAG_INFO << "Subtraction:   a - b = " << m << ", error = " << errorQVector(m, m_expc) << "\n";
    
    // Test the dot product:
    dcomplex d = a.dot(b);
    dcomplex d_expc = dcomplex(5.02392213375248,4.025019781801383);
    std::cout << TAG_INFO << "Dot product:   a·b   = " << d << ", error = " << abs(d - d_expc)/abs(d) << "\n";
    
    // Test the cross product:
    QVector c = a.cross(b);
    QVector c_expc = QVector(dcomplex(0.09823284143656252,7.763664533693396),
                             dcomplex(-5.502193196068165,-3.850646668302501),
                             dcomplex(3.873261330170582,0.4722029183180314));
    std::cout << TAG_INFO << "Cross product: a x b = " << c << ", error = " << errorQVector(c, c_expc) << "\n";
    
    return 0;
}

/**
 * Test the rotation of qvectors.
 */
int testRotations() {
    std::cout << "====== TEST ROTATIONS ======\n";
    
    // Test again the conversion from angle parameters to vector form:
    QVector a = QVector(dcomplex(0.7, 0.4), dcomplex(1.1, -0.5));
    QVector a_expc = QVector(dcomplex(0.7683955263452069,-0.4613632101591258),
                             dcomplex(-0.6964459431223345,-0.31416070729921625),
                             dcomplex(0.545810459257944,0.2486457969384876));
    std::cout << TAG_INFO << "a = " << a << ", error = " << errorQVector(a, a_expc) << "\n";
    
    // Test the rotation around the x axis:
    QVector arx = a.rotate_1(dcomplex(0.1, 0.5));
    QVector arx_expc = QVector(dcomplex(0.7683955263452069,-0.4613632101591258),
                               dcomplex(-0.7302742644199023,-0.6272445473990618),
                               dcomplex(0.7098180745615855,-0.1458840334837435));
    std::cout << TAG_INFO << "Rx.a = " << arx << ", error = " << errorQVector(arx, arx_expc) << "\n";
    std::cout << TAG_INFO << "norm(Rx.a) = " << arx.norm() << "\n";
    
    // Test the rotation around the y axis:
    QVector ary = a.rotate_2(dcomplex(-0.3, 1.2));
    QVector ary_expc = QVector(dcomplex(0.8843493039557755,0.1987394624213622),
                               dcomplex(-0.6964459431223345,-0.31416070729921625),
                               dcomplex(0.5790712545558938,-0.6813514783261888));
    std::cout << TAG_INFO << "Ry.a = " << ary << ", error = " << errorQVector(ary, ary_expc) << "\n";
    std::cout << TAG_INFO << "norm(Ry.a) = " << ary.norm() << "\n";
    
    // Test the rotation around the z axis:
    QVector arz = a.rotate_3(dcomplex(0.4, -0.8));
    QVector arz_expc = QVector(dcomplex(1.7258226657759863,-0.7086598449018777),
                               dcomplex(-0.7264704207204491,-1.4966992667906656),
                               dcomplex(0.545810459257944,0.2486457969384876));
    std::cout << TAG_INFO << "Rz.a = " << arz << ", error = " << errorQVector(arz, arz_expc) << "\n";
    std::cout << TAG_INFO << "norm(Rz.a) = " << arz.norm() << "\n";
    
    return 0;
}

/**
 * Test the shearing of qvectors (used in the contact interactions).
 */
int testShearing() {
    std::cout << "====== TEST SHEARING ======\n";
    
    QVector a = QVector(0.1 + 0.5*I, -0.5 + 0.3*I);
    
    QVector ashp = a.shear_plus(dcomplex(0.2, 0.5));
    QVector ashp_expc = QVector(dcomplex(-0.6795385786481603,-0.19034907140169013),
                                dcomplex(0.4036893212319249,-0.6496332202930164),
                                dcomplex(0.924227095441617,0.14379632124828884));
    std::cout << TAG_INFO << "Sp.a = " << ashp << ", error = " << errorQVector(ashp, ashp_expc) << "\n";
    std::cout << TAG_INFO << "norm(Sp.a) = " << ashp.norm() << "\n";
    
    QVector ashm = a.shear_minus(dcomplex(0.1, 0.6));
    QVector ashm_expc = QVector(dcomplex(-0.6747543511025345,1.1003991415396572),
                                dcomplex(0.6619093860969617,-0.3921350061988173),
                                dcomplex(1.4066610810235156,0.7123655892963552));
    std::cout << TAG_INFO << "Sp.a = " << ashm << ", error = " << errorQVector(ashm, ashm_expc) << "\n";
    std::cout << TAG_INFO << "norm(Sp.a) = " << ashm.norm() << "\n";
    
    return 0;
}

/**
 * Main function of the test of the qvector object.
 */
int main(int argc, char** argv) {
    
    basicOperations();
    testRotations();
    testShearing();
    
    return 0;
}
