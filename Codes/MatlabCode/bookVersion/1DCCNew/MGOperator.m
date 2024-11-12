function u = MGOperator(u,f,hf,level,MGParam)

  nfPlusGhostLayers = length(u);
  nf = nfPlusGhostLayers-2;

  u = smooth(u,f,hf,MGParam.m1,MGParam.omega);

  if level > 0

    hc = 2*hf;
    nc = nf/2;

    cGr = restriction(getResidual(u,f,hf));
    cGc = zeros(nc+2,1);

    for sig = 1:MGParam.pCycle
      cGc = MGOperator(cGc,cGr,hc,level-1,MGParam);
    end

    u(2:nf+1) = u(2:nf+1)+prolongation(cGc(2:nc+1));

    u = smooth(u,f,hf,MGParam.m2,MGParam.omega);
  end
  
end