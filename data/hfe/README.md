# HFE Data

Common folder to potentially hold data downloaded from HFE database.

The documentation ("HFE_ProductSpecifications.pdf") was obtained from
https://open.canada.ca/data/en/dataset/fe83a604-aa5a-4e46-903c-685f8b0cc33c.

The database in JSON format was provided by Philippe Aussant
(philippe.aussant@NRCan-RNCan.gc.ca) on Aug 29, 2022.

The following data were provided:
* ``historical_flood.json`` ... This layer contains individual flood
  occurrences (``points``). The points in this layer have a unique ID
  (``uuid``) and a secondary key to their historical_flood_events
  (``event_id``). We need precipitation for these points. You can
  write the precipitation for each in the ``rainfall_mm``
  attribute. The maps produced in 2.4.1 of SOW are for these points.
* ``historical_flood_event.json`` ... This layer contains the same
  flood occurrences as the previous layer but now they are grouped
  into events (``multipoints``). A flood event contains one or more
  occurrence. The unique id of this layer is ``event_id``. The maps
  produced in 2.4.2 and 2.4.3 of SOW are for these multipoints.
* ``large_historical_flood_event_sample.json`` ... As discussed in the
  kickoff meeting, this file contains a sample of the
  ``historical_flood_event`` layer. It has three large flood events
  that contain multiple flood occurrence.
* ``small_historical_flood_event_sample.json`` ... As discussed in the
  kickoff meeting, this file contains a sample of the
  ``historical_flood_event`` layer. It has three small flood events
  that contain only one occurrence.
