using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class entrance : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        var prefab = Resources.Load("AVE");
        var go = UnityEngine.Object.Instantiate(prefab, UnityEngine.Vector3.zero, Quaternion.identity) as GameObject;
        if (go == null)
        {
            Debug.Log("why ??????????");
        }
        else
        {
            Debug.Log(go.name);
            var defaultGo = go.transform.FindChild("default").gameObject;
            var comp = defaultGo.AddComponent<RotateModel>();
            comp.AVE = defaultGo;
            Debug.Log(comp.AVE == null ? "AVE NULL !!!!" : comp.AVE.name);
        }

    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
