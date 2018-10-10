    function [] = plot_turn_in_place(k, H, sgn, gamma, x, y, angle, s)
        c = [0.0, 0.0, 0.0];
        r = k*1.2;
        fih = H*sgn;

        c(1) = x;
        c(2) = y;
        c(3) = angle;
        fimax = c(3) + s*sgn;
           
        xx2(1) = c(1);
        xx2(2) = c(2);
        yy2(1) = c(1) + cos(gamma)*k;
        yy2(2) = c(2) + sin(gamma)*k;
        arrow(xx2, yy2, 'Length', 8, 'BaseAngle', 10, 'TipAngle', 20);
     
        fi = c(3) : fih : fimax;
        x3 = r * cos(fi) + c(1);
        y3 = r * sin(fi) + c(2);
        plot(x3, y3, 'r');  
        
        sz = size(x3);
        if (sz(2) >= 10)
            xx0(1) = x3(6);
            xx0(2) = y3(6);
            yy0(1) = x3(9);
            yy0(2) = y3(9);
            arrow(xx0, yy0, 'Length', 8, 'BaseAngle', 20, 'TipAngle', 10, 'EdgeColor','b','FaceColor','b');
        end;
    end