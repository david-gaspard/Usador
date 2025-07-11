/****
 * @date Created on 2025-07-03 at 14:35:12 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the methods for the Contact interaction object.
 ***/
#include "Contact.hpp"

/**
 * Constructor of the contact interaction.
 */
Contact::Contact(const Vector2D c0, const Vector2D c1, const int stype, const dcomplex ga) {
    if (stype != +1 && stype != -1) {
        throw std::invalid_argument("In Contact(): Invalid ltype, expected +1 or -1.");
    }
    this->c0 = c0;
    this->c1 = c1;
    this->stype = stype;
    this->ga = ga;
}

/**
 * Returns the value of gamma of this contact interaction.
 */
dcomplex Contact::getGamma() const {
    return ga;
}

/**
 * Assigns the value of gamma of this contact interaction.
 */
void Contact::setGamma(const dcomplex ga) {
    this->ga = ga;
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
