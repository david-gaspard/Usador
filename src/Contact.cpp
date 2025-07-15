/****
 * @date Created on 2025-07-03 at 14:35:12 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the methods for the Contact interaction object.
 ***/
#include "Contact.hpp"
#include "UsadelSystem.hpp"

/**
 * Constructor of the contact interaction.
 */
Contact::Contact(const Vector2D c0, const Vector2D c1, const int stype, const double tval) {
    if (stype != +1 && stype != -1) {
        throw std::invalid_argument("In Contact(): Invalid stype, expected +1 or -1.");
    }
    this->c0 = c0;
    this->c1 = c1;
    this->stype = stype;
    setTransmission(tval);
}

/**
 * Assigns the transmission value to the present Contact interaction object.
 */
void Contact::setTransmission(const double tval) {
    if (tval <= 0. || tval > 1.) {
        throw std::invalid_argument("In Contact(): Invalid transmission value, expected in (0, 1].");
    }
    this->tval = tval;
    ga = dcomplex(std::sqrt(1./tval), SQRTEPS);
}

/**
 * Returns the transmission value of the present Contact interaction object.
 */
double Contact::getTransmission() const {
    return tval;
}

/**
 * Returns the value of gamma of this contact interaction.
 */
dcomplex Contact::getGamma() const {
    return ga;
}

/**
 * Returns the shearing type of the contact interaction: +1 for Sigma_+, and -1 for Sigma_-.
 */
int Contact::getSType() const {
    return stype;
}

/**
 * Returns True when the given line segment p0-p1 intersects the Contact interaction.
 */
bool Contact::intersects(Vector2D p0, Vector2D p1) const {
    
    double det = (c1 - c0).cross(p1 - p0);
    
    if (det == 0.) return false;  // If the segments are parallel, then stops and returns False.
    
    double s = (p0 - c0).cross(p1 - p0)/det;  // Parameter along the contact inter. segment: c(s) = c0 + (c1-c0)*s
    double t = (p0 - c0).cross(c1 - c0)/det;  // Parameter along the propagation segment:    p(t) = p0 + (p1-p0)*t
    
    return (0. < s && s < 1. && 0. < t && t < 1.);
}

/**
 * Main computational routine of the contact interaction.
 * Applies a transformation to the QVector "qv" attached to point "p1" given that we propagate from point "p1" to the central point "p0".
 * If the contact is not crossed by p1-p0, then "qv" is left unchanged.
 * In general, "p0" is the central point of a five-point stencil, and "p1" is one of the four neighboring points.
 * 
 * See also: UsadelSystem::residualPoint()
 */
void Contact::applyTransfo(Vector2D p0, Vector2D p1, QVector& qv) const {
    
    double det = (c1 - c0).cross(p1 - p0);
    
    if (det != 0.) {// If the segments are not parallel, then continues.
    
        double s = (p0 - c0).cross(p1 - p0)/det;  // Parameter along the contact inter. segment: c(s) = c0 + (c1-c0)*s
        double t = (p0 - c0).cross(c1 - c0)/det;  // Parameter along the propagation segment:    p(t) = p0 + (p1-p0)*t
        
        if (0. < s && s < 1. && 0. < t && t < 1.) {// If "s" and "t" are in the unit interval [0, 1], then there is an intersection.
            
            double sign = det < 0. ? 1. : -1.;
            // If det < 0, then apply direct transformation: Q -> exp(i*gamma*Sigma_+-)*Q*exp(-i*gamma*Sigma_+-).
            // If det > 0, then apply opposite transformation.
            
            if (stype == +1) {// If the transformation involves Lambda_+ (Sigma_+).
                qv = qv.shear_plus(sign * I * ga);
            }
            else {// Otherwise, if the transformation involves Lambda_- (Sigma_-).
                qv = qv.shear_minus(sign * I * ga);
            }
            
        }
    }
}

/**
 * Returns the length of the contact interaction. This is the 2D equivalent of the surface area in 3D.
 */
double Contact::getLength() const {
    return (c1 - c0).norm();
}

/**
 * Returns the normal vector of the contact interaction in the direction of the flux, i.e., in the direction to which the direct transformation
 * Q -> exp(i*gamma*Sigma_+-)*Q*exp(-i*gamma*Sigma_+-) is applied.
 */
Vector2D Contact::getNormal() const {
    return Vector2D(c0.y - c1.y, c1.x - c0.x).normalize();
}

/**
 * Returns the flux integral: Integral tr( J.Lambda_+- ) dS , where J is the matrix current density J = -(lscat/d) Q Grad(Q).
 * 
 * @deprecated This is an old implementation unable to take into account the boundary conditions.
 */
dcomplex Contact::getFlux_old(const UsadelSystem& usys) const {
    
    const QVector lv(dcomplex(1., 0.), dcomplex(0., stype), dcomplex(0., 0.));  // Vector representation of Lambda_+-.
    
    MeshPoint p0, p1;
    QVector jv, flux;
    Vector2D dir, normal;
    
    normal = getNormal();  // Compute the normal vector to the contact interaction.
    std::cout << TAG_INFO << "Normal vector = " << normal << "\n";
    int nintersec = 0; // Initialize the number of intersections.
    
    // TODO: REWRITE THIS FUNCTION COMPLETELY to deal with contact interactions out of the mesh.
    // To this end, we have to check for neighboring points p0 + (1, 0), p0 + (-1, 0), p0 + (0, 1), p0 + (0, -1)...
    // This is needed to check the bimodal distribution [TEST: testRhoBimodal()]....
    // .............................................................................
    
    
    for (int i = 0; i < usys.getNPoint(); i++) {// Loop on the points.
        
        p0 = usys.getPoint(i);  // Get the point.
        
        if (p0.south >= 0) {// If the south point belongs to the mesh.
            p1 = usys.getPoint(p0.south);
            if (intersects(Vector2D(p0), Vector2D(p1))) {// If the segment p0-p1 intersect the contact interaction, then add the contribution to the flux.
                //std::cout << TAG_INFO << "Found intersection between p0 = " << Vector2D(p0) << " and p1 = " << Vector2D(p1) << " (south point), ";
                usys.getJVector(i, p0.south, jv, dir);
                //std::cout << "jv = " << jv << ", dir = " << dir << "\n";
                flux += jv * dir.dot(normal);
                nintersec++;
            }
        }
        if (p0.east >= 0) {// If the east point belongs to the mesh.
            p1 = usys.getPoint(p0.east);
            if (intersects(Vector2D(p0), Vector2D(p1))) {// If the segment p0-p1 intersect the contact interaction, then add the contribution to the flux.
                //std::cout << TAG_INFO << "Found intersection between p0 = " << Vector2D(p0) << " and p1 = " << Vector2D(p1) << " (east point), ";
                usys.getJVector(i, p0.east, jv, dir);
                //std::cout << "jv = " << jv << ", dir = " << dir << "\n";
                flux += jv * dir.dot(normal);
                nintersec++;
            }
        }
    }
    
    std::cout << TAG_INFO << "flux = " << flux << ", length = " << getLength() << ", nintersec = " << nintersec << "\n";
    
    flux *= getLength()/nintersec;  // Normalize the flux integral.
    
    return flux.dot(lv); // Returns the component tr(flux.Lambda_+-), instead of the matrix flux.
}


/**
 * Returns the flux integral: Integral tr( J.Lambda_+- ) dS , where J is the matrix current density J = -(lscat/d) Q Grad(Q).
 */
dcomplex Contact::getFlux(const UsadelSystem& usys) const {
    
    const QVector lv(dcomplex(1., 0.), dcomplex(0., stype), dcomplex(0., 0.));  // Vector representation of Lambda_+-.
    QVector jv, flux;
    MeshPoint p0;
    Vector2D p0v, dir, normal;
    int nintersec = 0; // Initialize the number of intersections.
    
    for (int i = 0; i < usys.getNPoint(); i++) {// Loop on the points.
        
        p0 = usys.getPoint(i);  // Get the point.
        p0v = Vector2D(p0);
        
        if (intersects(p0v, p0v + Vector2D(0., 1.))) {
            usys.getJVector(i, NORTH, jv, dir);  // NB: This assumes another implementation of the getJVector() method.
            flux += jv * normal.dot(Vector2D(0., 1.));
            nintersec++;
        }
        
        // WARNING: This loop crosses the contact TWO TIMES per segment !!!
        // TWO OPTIONS: 
        // 1. Loop over the segments uniquely using the previous technique (south, east), but taking into account the (north, west) directions at the edges of the medium.
        // 2. Loop directly over the crossing segments by looking for the ranges [xmin, xmax] and [ymin, ymax] containing the contact interaction.
        
    }
    
    return flux.dot(lv); // Returns the component tr(flux.Lambda_+-), instead of the matrix flux.
}
