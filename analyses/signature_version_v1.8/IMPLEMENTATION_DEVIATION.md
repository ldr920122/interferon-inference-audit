# Historical implementation note

The prestated Chinese note includes a constant background column. The archived numpy inference and official Util.ridge_significance_test calls did not set flag_const=True; official default is False at commit ce53cc7f. Thus the actual analyzed outputs do not include a constant column. The manuscript now states this discrepancy. The historical note is retained without retrospective alteration. No new intercept variant or tuned analysis was substituted into the primary results.
