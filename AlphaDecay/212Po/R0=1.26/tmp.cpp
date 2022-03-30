//    for (int i=1; i<N-1; i++){
//        /* N in paramsPTR is number of grid points. 
//         * Begin here at second gridpoint (i.e. R) and go to
//         * next to last gridpoint.*/
//         
//    
//        V = paramsPTR -> V(i);
//        r = paramsPTR -> r(i);
//
//        if (V > E){
//            // Expontetnial functions
//
//            omega = sqrt(alpha*(V-E));
//            // Continuity of function to the left
//            A(2*i-2,2*i-1) = -exp(omega*r);
//            A(2*i-2,2*i-0) = -exp(-omega*r);
//
//            // Continuity of deriv. to the left
//            A(2*i-1,2*i-1) = -omega*exp(omega*r);
//            A(2*i-1,2*i+0) = +omega*exp(-omega*r);
//
//            r = paramsPTR -> r(i+1);
//            // Continuity of function to the right
//            A(2*i+0,2*i-1) = +exp(omega*r);
//            A(2*i+0,2*i-0) = +exp(-omega*r);
//
//            // Continuity of deriv. to the right
//            A(2*i+1,2*i-1) = +omega*exp(omega*r);
//            A(2*i+1,2*i+0) = -omega*exp(-omega*r);
//            cout << '+';
//        } else { 
//
//        if (V < E) {
//            // Trigonometric functions
//
//            omega = sqrt(alpha*(E-V));
//            // Continuity of function to the left
//            A(2*i-2,2*i-1) = -cos(omega*r);
//            A(2*i-2,2*i-0) = -sin(-omega*r);
//
//            // Continuity of deriv. to the left
//            A(2*i-1,2*i-1) = +omega*sin(omega*r);
//            A(2*i-1,2*i+0) = -omega*cos(-omega*r);
//
//            r = paramsPTR -> r(i+1);
//            // Continuity of function to the right
//            A(2*i+0,2*i-1) = +cos(omega*r);
//            A(2*i+0,2*i-0) = +cos(-omega*r);
//
//            // Continuity of deriv. to the right
//            A(2*i+1,2*i-1) = -omega*sin(omega*r);
//            A(2*i+1,2*i+0) = +omega*cos(-omega*r);
//            cout << '-';
//        } 
//        else { 
//            /* Linear function. Haven't taken this
//             * into account. */ 
//            cout << "E=V at some point. Terminating..." << endl;    
//            exit(1);
//        }
//        }
//    }
//
//
//    V = paramsPTR -> V(N-1);
//    r = paramsPTR -> r(N-1);
//    if (V > E) {
//        // Exponential function
//
//        omega = sqrt(alpha*(V-E));
//        // Continuity of function to the left
//        A(2*(N-1)-2,2*(N-1)-1) = -exp(-omega*r);
//
//        // Continuity of deriv. to the left
//        A(2*(N-1)-1,2*(N-1)-1) = +omega*exp(-omega*r);
//    } else {
//    
//    if (V < E) {
//        // ~exp(i*omega*x)
//
//        omega = sqrt(alpha*(E-V));
//        
//        // Continuity of function to the left
//        A(2*(N-1)-2,2*(N-1)-1) = -(cos(omega*r),sin(omega*r));
//
//        // Continuity of deriv. to the left
//        A(2*(N-1)-1,2*(N-1)-1) = -omega*(-sin(omega*r),cos(omega*r));
