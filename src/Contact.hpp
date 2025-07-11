/****
 * @date Created on 2025-07-03 at 14:21:23 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the definition of the contact interaction object.
 ***/
#ifndef _CONTACT_H
#define _CONTACT_H
#include "Vector2D.hpp"
#include "QVector.hpp"

/**
 * Define the contact interaction.
 */
class Contact {
    
    private:
    
    Vector2D c0;  // Position of the first end of the contact interaction (usually above).
    Vector2D c1;  // Position of the second end of the contact interaction (usually below).
                  // Note that the normal vector (forward direction of the transformation) is (-dy, dx) where d = c1 - c0.
    int stype;    // Type of shearing of the contact interaction: +1 for Sigma_+, and -1 for Sigma_-.
    dcomplex ga;  // Gamma parameter of the contact interaction.
    
    public:
    
    // Constructors:
    Contact(const Vector2D c0, const Vector2D c1, int stype, const dcomplex ga);
    
    // Getters:
    dcomplex getGamma() const;
    
    // Setters:
    void setGamma(const dcomplex ga);
    
    // Computational methods:
    bool intersects(Vector2D p0, Vector2D p1) const;
    void applyTransfo(Vector2D p0, Vector2D p1, QVector& qv) const;
    
};

#endif
