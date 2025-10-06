function ParIndexes = SelectParents_Fcn(Cost,SelectionNum,SelMethod)
PopSize = size(Cost,1);

switch SelMethod
    case 1
        R = randperm(PopSize); ParIndexes = R(1:SelectionNum);
end

end