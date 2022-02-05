//
//  Rosembrok_function.h
//  Rosembrok_function
//
//  Created by Cibrán López Álvarez on 02/02/2022.
//  Definition of the Rosembrok function.


#ifndef Rosembrok_function_h
#define Rosembrok_function_h

#define a 1
#define b 100


double Rosembrok(double x, double y) {
    return b * (y - x * x) * (y - x * x) + (a - x) * (a - x);
}

double Rosembrok_dx(double x, double y) {
    /* Defines the derivative with respect to 'x' of the Rosembrok function */
    return - 4 * b * x * (y - x * x) - 2 * (a - x);
}

double Rosembrok_dy(double x, double y) {
    /* Defines the derivative with respect to 'y' of the Rosembrok function */
    return 2 * b * (y - x * x);
}

#endif /* Rosembrok_function_h */
