%% Allocate PNum size matrices here
variables_initializer=zeros(PNum,1);

%% initialize intializers
disp('initializing initializers')
for f=1:FNum
    x=tilde_node_position(f,1);
    y=tilde_node_position(f,2);
    variables_initializer(first_p_for_f(f))=A*exp(-((x-gauss_center_x)^2+(y-gauss_center_y)^2)/sigma) ...
        *b_area(f);
end