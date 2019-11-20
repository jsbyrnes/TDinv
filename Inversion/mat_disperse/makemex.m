%make the mex files

codegen modrt -args {complex(zeros(2,2,101)), complex(zeros(2,2,101)), complex(zeros(2,2,101)), complex(zeros(2,2,101)), complex(zeros(2,2,100))}
codegen genrt -args {complex(zeros(2,2,100)),complex(zeros(2,2,100)),complex(zeros(2,2,100)),complex(zeros(2,2,100))}
codegen psv -args { zeros(100,1), zeros(101,1), zeros(101,1), zeros(101,1), 0,0 }

% codegen modrt -args {complex(zeros(2,2,81)), complex(zeros(2,2,81)), complex(zeros(2,2,81)), complex(zeros(2,2,81)), complex(zeros(2,2,80))}
% codegen genrt -args {complex(zeros(2,2,80)),complex(zeros(2,2,80)),complex(zeros(2,2,80)),complex(zeros(2,2,80))}
% codegen psv -args { zeros(80,1), zeros(81,1), zeros(81,1), zeros(81,1), 0,0 }