//
//  Deterministic_optimization.c
//  Rosembrok_function
//
//  Created by Cibrán López Álvarez on 20/01/2022.
//  Performance of the steepest descent and conjugate gradient methods.


#include <stdio.h>
#include "Rosembrok_function.h"

#define x_0 -1.5
#define y_0 -1
#define tolerance_square 1e-14 // Square, so no unnecessary square roots are calculated

#define maximum_iterations 1000
#define rho_backtracking 0.996
#define sigma_backtracking 0.501
#define n_partial 10


int not_converged(double D_x, double D_y, unsigned int iteration) {
    /* The optimization will be performed until the modulus of the gradient is lower than tolerance or a maximum number of steps is reached */
    if ((D_x*D_x + D_y*D_y) <= tolerance_square)
        return 0;
    
    if (iteration >= maximum_iterations)
        return 0;
    
    return 1;
}

void backtracking_line_search(double *alpha, double (*f)(double, double), double x, double y, double d_x , double d_y, double D_x , double D_y) {
    /* Obtains alpha such that the new minimum point... */
    *alpha = 1;
    
    while ((*f)(x + *alpha * D_x, y + *alpha * D_y) > (*f)(x, y) + sigma_backtracking * *alpha * (d_x * D_x + d_y * D_y))
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

void steepest_descent_method(double (*f)(double, double), double (*f_dx)(double, double), double (*f_dy)(double, double)) {
    unsigned int iteration = 0;
    double x, y, d_x, d_y, D_x, D_y, alpha;
        
    x = x_0; y = y_0;
    d_x = (*f_dx)(x, y); d_y = (*f_dy)(x, y);
    D_x = - d_x; D_y = - d_y;
    
    while (not_converged(d_x, d_y, iteration)) {
        iteration++;
        
        /* Finding a good parameter for the next step */
        backtracking_line_search(&alpha, (*f), x, y, d_x, d_y, D_x, D_y);
        
        /* Advancing one step forward */
        x += alpha * D_x;
        y += alpha * D_y;
        
        d_x = (*f_dx)(x, y);
        d_y = (*f_dy)(x, y);
        
        /* Updating the descending direction */
        D_x = - d_x;
        D_y = - d_y;
        printf("%f\n", Rosembrok(-1.5, -1));
    }
    
    printf("With the steepest descent method, it converges to (x, y) = (%.15lf, %.15lf) in %u iterations, with f(x, y) = %.15lf\n", x, y, iteration, (*f)(x, y));
}

void conjugate_gradient_method(double (*f)(double, double), double (*f_dx)(double, double), double (*f_dy)(double, double), void (*beta_f)(double*, double, double, double, double), int partial) {
    unsigned int iteration = 0;
    double x, y, d_x, d_y, D_x, D_y, alpha, beta = 0;
    
    x = x_0; y = y_0;
    d_x = (*f_dx)(x, y); d_y = (*f_dy)(x, y);
    D_x = - d_x; D_y = - d_y;
    
    while (not_converged(d_x, d_y, iteration)) {
        iteration++;
        
        /* Finding a good alpha parameter for the next step */
        backtracking_line_search(&alpha, (*f), x, y, d_x, d_y, D_x, D_y);
        
        /* Advancing one step forward */
        x += alpha * D_x;
        y += alpha * D_y;
        
        /* Calculating beta parameter */
        (*beta_f)(&beta, (*f_dx)(x, y), (*f_dy)(x, y), d_x, d_y);
        
        if (partial && ((iteration % n_partial) == 0))
            beta = 0;
        //printf("%f", beta);
        /* Saving the previous derivative */
        d_x = (*f_dx)(x, y);
        d_y = (*f_dy)(x, y);
        
        /* Updating the descending direction */
        D_x = beta * D_x - d_x;
        D_y = beta * D_y - d_y;
    }
    
    printf("With the ");
    
    if (partial)
        printf("partial ");
    
    printf("conjugate gradient method, it converges to (x, y) = (%.15lf, %.15lf) in %u iterations, with f(x, y) = %.15lf\n", x, y, iteration, (*f)(x, y));
}


int main() {
    /* Here we call the main funtions, which have all the neccessary parameters defined (so it is easier to work separataly with them) */
    
    steepest_descent_method(Rosembrok, Rosembrok_dx, Rosembrok_dy);
    
    printf("\nFletcher Reeves method:\n");
    
    conjugate_gradient_method(Rosembrok, Rosembrok_dx, Rosembrok_dy, Fletcher_Reeves_method, 0);
    
    conjugate_gradient_method(Rosembrok, Rosembrok_dx, Rosembrok_dy, Fletcher_Reeves_method, 1);
    
    printf("\nPolak Ribiere method:\n");
    
    conjugate_gradient_method(Rosembrok, Rosembrok_dx, Rosembrok_dy, Polak_Ribiere_method, 0);
    
    conjugate_gradient_method(Rosembrok, Rosembrok_dx, Rosembrok_dy, Polak_Ribiere_method, 1);
}
//
//  rosembrok_function.c
//  Rosembrok_function
//
//  Created by Cibrán López Álvarez on 20/01/2022.
//  Performance of the steepest descent and conjugate gradient methods.


#include <stdio.h>
#include "Rosembrok_function.h"

#define x_0 -1.5
#define y_0 -1
#define tolerance_square 1e-14 // Square, so no unnecessary square roots are calculated

#define maximum_iterations 1000
#define rho_backtracking 0.996
#define sigma_backtracking 0.501
#define n_partial 10


int not_converged(double D_x, double D_y, unsigned int iteration) {
    /* The optimization will be performed until the modulus of the gradient is lower than tolerance or a maximum number of steps is reached */
    if ((D_x*D_x + D_y*D_y) <= tolerance_square)
        return 0;
    
    if (iteration >= maximum_iterations)
        return 0;
    
    return 1;
}

void backtracking_line_search(double *alpha, double (*f)(double, double), double x, double y, double d_x , double d_y, double D_x , double D_y) {
    /* Obtains alpha such that the new minimum point... */
    *alpha = 1;
    
    while ((*f)(x + *alpha * D_x, y + *alpha * D_y) > (*f)(x, y) + sigma_backtracking * *alpha * (d_x * D_x + d_y * D_y))
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

void steepest_descent_method(double (*f)(double, double), double (*f_dx)(double, double), double (*f_dy)(double, double)) {
    unsigned int iteration = 0;
    double x, y, d_x, d_y, D_x, D_y, alpha;
        
    x = x_0; y = y_0;
    d_x = (*f_dx)(x, y); d_y = (*f_dy)(x, y);
    D_x = - d_x; D_y = - d_y;
    
    while (not_converged(d_x, d_y, iteration)) {
        iteration++;
        
        /* Finding a good parameter for the next step */
        backtracking_line_search(&alpha, (*f), x, y, d_x, d_y, D_x, D_y);
        
        /* Advancing one step forward */
        x += alpha * D_x;
        y += alpha * D_y;
        
        d_x = (*f_dx)(x, y);
        d_y = (*f_dy)(x, y);
        
        /* Updating the descending direction */
        D_x = - d_x;
        D_y = - d_y;
        printf("%f\n", Rosembrok(-1.5, -1));
    }
    
    printf("With the steepest descent method, it converges to (x, y) = (%.15lf, %.15lf) in %u iterations, with f(x, y) = %.15lf\n", x, y, iteration, (*f)(x, y));
}

void conjugate_gradient_method(double (*f)(double, double), double (*f_dx)(double, double), double (*f_dy)(double, double), void (*beta_f)(double*, double, double, double, double), int partial) {
    unsigned int iteration = 0;
    double x, y, d_x, d_y, D_x, D_y, alpha, beta = 0;
    
    x = x_0; y = y_0;
    d_x = (*f_dx)(x, y); d_y = (*f_dy)(x, y);
    D_x = - d_x; D_y = - d_y;
    
    while (not_converged(d_x, d_y, iteration)) {
        iteration++;
        
        /* Finding a good alpha parameter for the next step */
        backtracking_line_search(&alpha, (*f), x, y, d_x, d_y, D_x, D_y);
        
        /* Advancing one step forward */
        x += alpha * D_x;
        y += alpha * D_y;
        
        /* Calculating beta parameter */
        (*beta_f)(&beta, (*f_dx)(x, y), (*f_dy)(x, y), d_x, d_y);
        
        if (partial && ((iteration % n_partial) == 0))
            beta = 0;
        //printf("%f", beta);
        /* Saving the previous derivative */
        d_x = (*f_dx)(x, y);
        d_y = (*f_dy)(x, y);
        
        /* Updating the descending direction */
        D_x = beta * D_x - d_x;
        D_y = beta * D_y - d_y;
    }
    
    printf("With the ");
    
    if (partial)
        printf("partial ");
    
    printf("conjugate gradient method, it converges to (x, y) = (%.15lf, %.15lf) in %u iterations, with f(x, y) = %.15lf\n", x, y, iteration, (*f)(x, y));
}


int main() {
    /* Here we call the main funtions, which have all the neccessary parameters defined (so it is easier to work separataly with them) */
    
    steepest_descent_method(Rosembrok, Rosembrok_dx, Rosembrok_dy);
    
    printf("\nFletcher Reeves method:\n");
    
    conjugate_gradient_method(Rosembrok, Rosembrok_dx, Rosembrok_dy, Fletcher_Reeves_method, 0);
    
    conjugate_gradient_method(Rosembrok, Rosembrok_dx, Rosembrok_dy, Fletcher_Reeves_method, 1);
    
    printf("\nPolak Ribiere method:\n");
    
    conjugate_gradient_method(Rosembrok, Rosembrok_dx, Rosembrok_dy, Polak_Ribiere_method, 0);
    
    conjugate_gradient_method(Rosembrok, Rosembrok_dx, Rosembrok_dy, Polak_Ribiere_method, 1);
}
