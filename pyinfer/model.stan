functions{

    //*********************************************************************
    // Floor that returns an int
    //*********************************************************************    
  
    int floor_int(real x){
        int i;
        

        if (x>0){
            i=0;
            //while (i<x-1){i=i+1;}
            while (i<floor(x)){i=i+1;}
            
        }
        
        if (x<0){
            i=0;
            while (i>x){i=i-1;}
        }
        
        if (x==0){i=0;}
        
        return(i);
        
    
    }
    
    //*********************************************************************
    // Ceil that returns an int
    //*********************************************************************    
  
    int ceil_int(real x){
        int i;
        
        if (x>0){
            i=0;
            while (i<x){i=i+1;}
        }
        
        if (x<0){
            i=0;
            while (i>x+1){i=i-1;}
        }
        
        if (x==0){i=0;}
        
        return(i);
        
    
    }



    //*********************************************************************
    // 1D cubic interpolation
    //*********************************************************************

    real cubic_interpolator(array[] real p, real x){

        real val;

        val = p[1+1] + 0.5*x*(
                    p[2+1] - p[0+1] + x*(
                        2.0*p[0+1] - 5.0*p[1+1] + 4.0*p[2+1] - p[3+1] + x*(
                            3.0*(
                                p[1+1] - p[2+1]
                            ) + p[3+1] - p[0+1]
                        )
                    )
                );
                
        return val;
    }


    //*********************************************************************
    // 1D cubic 
    //*********************************************************************         
    
    real cubic(array[] real cube, real x, real x_min, real x_delta){
    
        array[4] real p;
        real r_x;
        int i_x;
        real d_x;
        real val;
        
        
        r_x = (x - x_min) / x_delta;
        
        i_x = floor_int( r_x  );
        
        d_x = r_x - floor(r_x);
        
        p = cube[ i_x-1+1 : i_x+2+1 ];  
        
        val = cubic_interpolator(p, d_x);
        
        return val;
    
    
    }

    //*********************************************************************
    // 1D cubic multi
    //*********************************************************************         
    

    
    array[] real cubic_multi(int N_chans, array[,] real cube, real x, real x_min, real x_delta){
    
        array[4] real p;
        real r_x;
        int i_x;
        real d_x;
        array[N_chans] real val;
        
        
        r_x = (x - x_min) / x_delta;
        
        i_x = floor_int( r_x  );
        
        d_x = r_x - floor(r_x);
        
        /*print(
            "x: ",x,
            "r_x: ",r_x,
            "i_x: ", i_x,
            "d_x: ", d_x);*/
        
        for (i in 1:N_chans){
    
            p = cube[i, i_x-1+1 : i_x+2+1 ];  

            val[i] = cubic_interpolator(p, d_x);

        }
        
        return val;
    
    
    }   

    

   //*********************************************************************
    // 1D cubic multi multi
    //*********************************************************************         
    

    
    //array[,] real cubic_all(int N_stn, int N_hrs, array[,,] real cube, real x, real x_min, real x_delta){
    matrix cubic_all(int N_stn, int N_hrs, array[,,] real cube, real x, real x_min, real x_delta){
    
        array[4] real p;
        real r_x;
        int i_x;
        real d_x;
        matrix[N_stn, N_hrs] val;
        
        
        r_x = (x - x_min) / x_delta;
        
        i_x = floor_int( r_x  );
        
        d_x = r_x - floor(r_x);
        
        for (i in 1:N_stn){
        
            for (j in 1:N_hrs){
            
                p = cube[i, j, i_x-1+1 : i_x+2+1 ];  
    
                val[i, j] = cubic_interpolator(p, d_x);

            }

        }
        
        return val;
    }

    //*********************************************************************
    // Linear interpolation from https://discourse.mc-stan.org/t/linear-interpolation-and-searchsorted-in-stan/13318/6
    //*********************************************************************         
    

    real linear_interpolation_v(real x, vector x_pred, vector y_pred){
        int K = rows(x_pred);
        vector[K] deltas = x - x_pred;
        real ans;
        real t;
        real w;
        int i;
        
        if(x<x_pred[1] || x>x_pred[K]) reject("x is outside of the x_pred grid!");
        
        if(rows(y_pred) != K) reject("x_pred and y_pred aren't of the same size");
        
        //this is which.min()
        //i = sort_indices_asc(fabs(deltas))[1];
    
        //
        i = sort_indices_asc(abs(deltas))[1];
    
        
        if(deltas[i]<=0) i -= 1;
        ans = y_pred[i];
        real x1 = x_pred[i];
        real x2 = x_pred[i + 1];
        real y1 = y_pred[i];
        real y2 = y_pred[i + 1];
        ans = y1 + (y2-y1) * (x-x1)/(x2-x1);
        t = (x-x1)/(x2-x1);
        
        w = 1-t;
        //w = 1/(1 + exp(1/(1-t) - 1/t));
        //w = 1 - 3*pow(t,2) + 2*pow(t,3);
        ans = w*y1 + (1-w)*y2;
        return ans;
  }

   //*********************************************************************
    // Linear multi
    //*********************************************************************         
    


    matrix linear_all_old(int N_stn, int N_hrs, array[,,] real cube, real x, real x_min, real x_delta){
    
        array[4] real p;
        real r_x;
        int i_x;
        real d_x;
        matrix[N_stn, N_hrs] val;
        real w;

       


        
        r_x = (x - x_min) / x_delta;
        
        i_x = floor_int( r_x  );
        
        d_x = r_x - floor(r_x);
        
        w = 1-d_x;
        //w = 1/(1 + exp(1/(1-d_x) - 1/d_x));
        //w = 1 - 3*pow(d_x,2) + 2*pow(d_x,3);

        
        for (i in 1:N_stn){
        
            for (j in 1:N_hrs){
            

                val[i, j] = cube[i,j,i_x+1]*(1-d_x) + cube[i,j,i_x+2]*d_x;

            }

        }
        
        return val;
    }




    matrix linear_all(int N_stn, int N_hrs, array[,,] real cube, real x, real x_min, real x_delta){
    
        array[4] real p;
        real r_x;
        int i_x;
        real d_x;
        matrix[N_stn, N_hrs] val;
        real w;
        
        r_x = (x - x_min) / x_delta;
        
        i_x = floor_int( r_x  );
        
        d_x = r_x - floor(r_x);
        
        w = 1-d_x;
        //w = 1/(1 + exp(1/(1-d_x) - 1/d_x));
        //w = 1 - 3*pow(d_x,2) + 2*pow(d_x,3);

        val = to_matrix(cube[:,:,i_x+1]) * (1-d_x) + to_matrix(cube[:,:,i_x+2]) * d_x;
        
        return val;
    }

    
}

data {



    int N_hrs;
    int N_stn;
    int N_puf;
    int N_hts;
    int N_oc;

    array[N_puf, N_stn, N_hrs, N_hts] real puffs;
    real height_min;
    real height_delta;

    
    vector[N_oc] obs_conc;
    array[N_oc, 2] int ij_conc;

    
    real height_lower;
    real height_upper;
        
    real flux_lower;
    real flux_upper;
    
    real sigma_conc_lower;
    real sigma_conc_upper;

    real sigma_conc_loc;
    real sigma_conc_scale;

    real flux_sigma_loc;
    real flux_sigma_scale;

    real flux_mu_loc;
    real flux_mu_scale;

    real height_sigma_loc; 
    real height_sigma_scale;

    real height_mu_loc;
    real height_mu_scale;


    int LIKELIHOOD;


}



parameters {

//https://mc-stan.org/docs/2_19/reference-manual/univariate-data-types-and-variable-declarations.html


    real<
        offset = height_mu_loc,
        multiplier = height_mu_scale
        > height_mu;

    real<
        offset = height_sigma_loc,
        multiplier = height_sigma_scale
        > height_sigma;

    real<
        offset = flux_mu_loc,
        multiplier = flux_mu_scale
        > flux_mu;

    real<
        offset = flux_sigma_loc,
        multiplier = flux_sigma_scale
        > flux_sigma;

    vector<
        lower = (height_lower - height_mu)/height_sigma, 
        upper = (height_upper - height_mu)/height_sigma
        >[N_puf] height0;
        
    vector<
        lower = (flux_lower - flux_mu)/flux_sigma, 
        upper = (flux_upper - flux_mu)/flux_sigma
        >[N_puf] flux0;

        
    real<
        lower = (sigma_conc_lower - sigma_conc_loc)/sigma_conc_scale, 
        upper = (sigma_conc_upper - sigma_conc_loc)/sigma_conc_scale
        > sigma_conc0;

}   

transformed parameters {

    vector[N_puf] flux;
    vector[N_puf] height;
    real sigma_conc;
    
    matrix[N_stn, N_hrs] puffs_interpolated;
    matrix[N_stn, N_hrs] conc_matrix;
    vector[N_oc] conc_vector;

    flux = flux0*flux_sigma + flux_mu;

    height = height0*height_sigma + height_mu;

    sigma_conc = sigma_conc0 * sigma_conc_scale + sigma_conc_loc;

    

    conc_matrix  = rep_matrix(0.0,N_stn, N_hrs); 


    for (i in 1:N_puf ){

            puffs_interpolated = linear_all(N_stn, N_hrs, puffs[i], height[i], height_min, height_delta);
            
            conc_matrix = conc_matrix + flux[i] * puffs_interpolated;

        
        }


    

    for (ii in 1:N_oc){
        int i;
        int j;

        i = ij_conc[ii,1];
        j = ij_conc[ii,2];
    
        conc_vector[ii] = conc_matrix[i,j];
    }


    

}

model {

    // PRIORS

    flux_sigma ~ normal(
                flux_sigma_loc,
                flux_sigma_scale
                );
    
    flux_mu ~ normal(
                flux_mu_loc,
                flux_mu_scale
                );
    
    height_sigma ~ normal(
                height_sigma_loc, 
                height_sigma_scale
                );

    height_mu ~ normal(
                height_mu_loc, 
                height_mu_scale
                );
    
    height0 ~ std_normal();

    flux0 ~ std_normal();

    sigma_conc0 ~ std_normal();

    // LIKELIHOOD

    if (LIKELIHOOD==1){    
        obs_conc ~ normal( conc_vector, sigma_conc);
    }
    


}

generated quantities{

    vector[N_oc] obs_hat;
    matrix[N_stn, N_hrs] conc_matrix_hat;


    for (i in 1:N_oc){

        obs_hat[i] = normal_rng( conc_vector[i], sigma_conc);
    
    }

    for (i in 1:N_stn){
        for (j in 1:N_hrs){
            conc_matrix_hat[i,j] = normal_rng( conc_matrix[i,j], sigma_conc);

        }
    
    }




}

