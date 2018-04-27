function [model] = model_EC_syms()

    %% Set parametrization
    model.param = 'log10';

    %% Set states
    syms SUBSTRATE ENZYME SUBSTRATE_ENZYME PRODUCT_ENZYME PRODUCT
    model.sym.x = [SUBSTRATE ENZYME SUBSTRATE_ENZYME PRODUCT_ENZYME PRODUCT];

    %% Set parameters
    syms enz_bind_fwd enz_bind_rev complex_trans_fwd complex_trans_rev prod_release_fwd prod_release_rev init_enzyme
    model.sym.p = [enz_bind_fwd enz_bind_rev complex_trans_fwd complex_trans_rev prod_release_fwd prod_release_rev init_enzyme];

    %% Set constants
    syms init_substrate
    model.sym.k = init_substrate;

    %% Set equations
    model.sym.xdot = sym(zeros(size(model.sym.x)));
    model.sym.xdot(1) = -enz_bind_fwd*SUBSTRATE*ENZYME + enz_bind_rev*SUBSTRATE_ENZYME;
    model.sym.xdot(2) = -enz_bind_fwd*SUBSTRATE*ENZYME + enz_bind_rev*SUBSTRATE_ENZYME + prod_release_fwd*PRODUCT_ENZYME - prod_release_rev*SUBSTRATE_ENZYME*PRODUCT;
    model.sym.xdot(3) =  enz_bind_fwd*SUBSTRATE*ENZYME - enz_bind_rev*SUBSTRATE_ENZYME + complex_trans_rev*PRODUCT_ENZYME - complex_trans_fwd*SUBSTRATE_ENZYME;
    model.sym.xdot(4) = -prod_release_fwd*PRODUCT_ENZYME + prod_release_rev*SUBSTRATE_ENZYME*PRODUCT - complex_trans_rev*PRODUCT_ENZYME + complex_trans_fwd*SUBSTRATE_ENZYME;
    model.sym.xdot(5) =  prod_release_fwd*PRODUCT_ENZYME - prod_release_rev*SUBSTRATE_ENZYME*PRODUCT;

    %% Set initial conditions
    model.sym.x0 = sym(zeros(size(model.sym.x)));
    model.sym.x0(1) = init_substrate;
    model.sym.x0(2) = init_enzyme;

    %% Set observables
    model.sym.y = sym(zeros(2,1));
    model.sym.y(1) = SUBSTRATE;
    model.sym.y(2) = PRODUCT;

end