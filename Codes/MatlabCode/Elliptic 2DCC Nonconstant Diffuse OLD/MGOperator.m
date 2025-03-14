function u = MGOperator(u,f,hf,level,MGParam)
  nfPlusGhostLayers = size(u);
  nf = nfPlusGhostLayers-2;

  u = smooth(u,f,hf,MGParam);

  if level > 0

    hc = 2*hf;
    nc = nf/2;

    cGr = restriction(getResidual(u,f,hf,MGParam));
    cGc = zeros(nc(1)+2,nc(2)+2);

    for sig = 1:MGParam.pCycle
      cGc = MGOperator(cGc,cGr,hc,level-1,MGParam);
    end

    u(2:nf(1)+1,2:nf(2)+1) = u(2:nf(1)+1,2:nf(2)+1) ...
      + prolongation(cGc(2:nc(1)+1,2:nc(2)+1));

    u = smooth(u,f,hf,MGParam);
  end
end