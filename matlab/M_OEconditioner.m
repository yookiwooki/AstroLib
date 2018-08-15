function conditioned_M_OE = M_OEconditioner(M_OE,mu)

for i = 1:length(M_OE)

COE(i,:) = M_OE2COE(M_OE(i,:));
RV(i,:) = COE2RV(COE(i,:),mu);
conditioned_COE(i,:) = RV2COE(RV(i,:));
conditioned_M_OE(i,:) = COE2M_OE(conditioned_COE(i,:));


end


end

