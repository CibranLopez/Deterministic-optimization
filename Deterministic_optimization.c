//
//  main.c
//  Rosembrok_function
//
//  Created by Cibrán López Álvarez on 20/01/2022.
//


#include <stdio.h>

#define x_0 -1.5
#define y_0 -1
#define tolerance 1e-14 // Square, so no unnecessary square roots are calculated

#define a 1
#define b 100

#define maximum_iterations 500
#define rho_backtracking 0.99
#define sigma_backtracking 0.5
#define n_partial 10



double Rosembrok(double x, double y) {
    /* Defines the Rosembrok function */
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

int not_converged(double D_x, double D_y, unsigned int iteration) {
    /* The optimization will be performed until the modulus of the gradient is lower than tolerance or a maximum number of steps is reached */
    if ((D_x*D_x + D_y*D_y) <= tolerance)
        return 0;
    
    if (iteration >= maximum_iterations)
        return 0;
    
    return 1;
}

void backtracking_line_search(double *alpha, double x, double y, double d_x , double d_y, double D_x , double D_y) {
    /* Obtains alpha such that the new minimum point... */
    *alpha = 1;
    
    while (Rosembrok(x + *alpha * D_x, y + *alpha * D_y) > Rosembrok(x, y) + sigma_backtracking * *alpha * (d_x * D_x + d_y * D_y))
        *alpha *= rho_backtracking;
}

void Fletcher_Reeves_method(double *beta, double d_x, double d_y, double old_dx, double old_dy) {
    /* Obtains beta according to the Fletcher Reeves method */
    *beta = (d_x * d_x + d_y * d_y) / (old_dx * old_dx + old_dy * old_dy);
}

void Polak_Ribiere_method(double *beta, double d_x, double d_y, double old_dx, double old_dy) {
    /* Obtains beta according to the Polak Ribiere method */
    *beta = ((d_x - old_dx) * d_x + (d_y - old_dy) * d_y) / (old_dx * old_dx + old_dy * old_dy);
}

void steepest_descent_method(void) {
    unsigned int iteration = 0;
    double x, y, d_x, d_y, D_x, D_y, alpha;
        
    x = x_0; y = y_0;
    d_x = Rosembrok_dx(x, y); d_y = Rosembrok_dy(x, y);
    D_x = - d_x; D_y = - d_y;
    
    while (not_converged(d_x, d_y, iteration)) {
        iteration++;
        
        /* Finding a good parameter for the next step */
        backtracking_line_search(&alpha, x, y, d_x, d_y, D_x, D_y);
        
        /* Advancing one step forward */
        x += alpha * D_x;
        y += alpha * D_y;
        
        d_x = Rosembrok_dx(x, y);
        d_y = Rosembrok_dy(x, y);
        
        /* Updating the descending direction */
        D_x = - d_x;
        D_y = - d_y;
        
        //printf("alpha: %f, iter: %u, x: %.15lf, y: %.15lf, f: %.15lf\n", alpha, iteration, x, y, Rosembrok(x, y));
    }
    
    printf("With the steepest descent method, it converges to (x, y) = (%.15lf, %.15lf) in %u iterations, with f(x, y) = %.15lf\n", x, y, iteration, Rosembrok(x, y));
}

void conjugate_gradient_method(int partial) {
    unsigned int iteration = 0;
    double x, y, d_x, d_y, D_x, D_y, alpha, beta = 0;
    
    x = x_0; y = y_0;
    d_x = Rosembrok_dx(x, y); d_y = Rosembrok_dy(x, y);
    D_x = - d_x; D_y = - d_y;
    
    while (not_converged(d_x, d_y, iteration)) {
        iteration++;
        
        /* Finding a good alpha parameter for the next step */
        backtracking_line_search(&alpha, x, y, d_x, d_y, D_x, D_y);
        
        /* Advancing one step forward */
        x += alpha * D_x;
        y += alpha * D_y;
        
        /* Calculating beta parameter */
        Fletcher_Reeves_method(&beta, Rosembrok_dx(x, y), Rosembrok_dy(x, y), d_x, d_y);
        
        if (partial && ((iteration % n_partial) == 0))
            beta = 0;
        
        /* Saving the previous derivative */
        d_x = Rosembrok_dx(x, y);
        d_y = Rosembrok_dy(x, y);
        
        /* Updating the descending direction */
        D_x = beta * D_x - d_x;
        D_y = beta * D_y - d_y;
        
        //printf("alpha: %f, beta: %f, iter: %u, x: %.15lf, y: %.15lf, f: %.15lf\n", alpha, beta, iteration, x, y, Rosembrok(x, y));
    }
    
    printf("With the ");
    
    if (partial)
        printf("partial ");
    
    printf("conjugate gradient method, it converges to (x, y) = (%.15lf, %.15lf) in %u iterations, with f(x, y) = %.15lf\n", x, y, iteration, Rosembrok(x, y));
}


int main() {
    /* Here we call the main funtions, which have all the neccessary parameters defined (so it is easier to work with them */
    
    steepest_descent_method();
    
    conjugate_gradient_method(0);
    
    conjugate_gradient_method(1);
}
