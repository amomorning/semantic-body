using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RotateModel : MonoBehaviour
{
    public GameObject AVE;
    Vector3[] vert;
    Mesh mesh;

    void Start()
    {
        AVE.gameObject.transform.Rotate(new Vector3(270, 315, 0));
        mesh = GetComponent<MeshFilter>().mesh;
        vert = mesh.vertices;
        float min = 0;
        for(var i = 0; i < vert.Length; ++ i)
        {
            float tmp = vert[i].z;
            if (tmp < min) min = tmp;
        }
        AVE.gameObject.transform.position += new Vector3(0.0f,-1.2f-min,0.0f);
    }

    void Update()
    {
        if (Input.GetMouseButton(0))
        {
            AVE.gameObject.transform.Rotate(new Vector3(0, 0, -Input.GetAxis("Mouse X") * Time.deltaTime * 300 ));
        }
    }
}
